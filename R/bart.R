# Creating the main BART function and it's generalities
new_tree <- function(x_train,
                     x_test){

        # Creating the list of node
        node_0 <- list(index = 0,
                       obs_train = 1:nrow(x_train),
                       obs_test  = 1:nrow(x_test),
                       left = NA,
                       right = NA,
                       parent = NA,
                       terminal = 1,
                       nog = 0,
                       depth = 0,
                       var = NA,
                       var_split_rule = NA,
                       mu = 0)
        tree <- list(node_0)
        names(tree) <- "node_0"
        # Returning the new tree
        return(tree)
}

# Getting terminal nodes
get_terminals <- function(tree){
       return( tree[unlist(lapply(tree, function(t){is.na(t$left)&is.na(t$right)}))] )
}

get_nonterminals <- function(tree){
        return( tree[unlist(lapply(tree, function(t){!is.na(t$left)&!is.na(t$right)}))] )
}

get_nog <- function(tree){
        return(tree[(unlist(lapply(tree,function(t){t$nog==1 & t$terminal==0 })))])
}

count_nog <- function(tree){
        return(sum(unlist(lapply(tree,function(t){t$nog==1 & t$terminal==0}))))
}

# Getting all nodes index
get_indexes <- function(tree){
        return(unlist(lapply(tree, function(t)(t$index))))
}

# Calculating the node likelihood
node_loglikelihood <- function(node,
                               res_vec,
                               tau,
                               tau_mu){

        # Slicing the current res_vec
        res_node <- res_vec[node$obs_train]
        n_obs <- length(res_node)

        return(-0.5*tau*crossprod(res_node)[[1]]-0.5*log(tau_mu+(n_obs*tau)) + (0.5*(tau^2)*(sum(res_node)^2))/(tau*n_obs+tau_mu))

}

# Creating a function to grow a tree
grow <- function(res_vec,
                 tree,
                 x_train,
                 x_test,
                 xcut,
                 tau,
                 tau_mu,
                 alpha,
                 beta,
                 node_min_size,
                 cat_var){

        # Getting the terminal nodes
        terminal_nodes <- get_terminals(tree)

        # Sampling one terminal node
        g_node_position <- sample(1:length(terminal_nodes),size = 1)
        g_node <- terminal_nodes[[g_node_position]]
        g_node_name <- names(terminal_nodes[g_node_position])
        g_node_position_orig <- which(names(tree) ==g_node_name)

        # Initializing the sample
        split_var_candidates <- colnames(x_train)
        good_tree_index <- 0

        # GETTING CATEGORICAL VARIABLES
        if(is.null(cat_var)){
                cat_var <- names(x_train)[!apply(x_train,2,function(x){length(unique(x))>round(nrow(x_train)*0.4)})]
        }

        while(good_tree_index==0){
                # Selecting a valid split
                split_var <- sample(split_var_candidates,size = 1)


                # For the case of categorical variables
                if(split_var %in% cat_var){

                        # Selecting the split rule for categorical variables
                        split_var_sampled_rule <- sample(unique(x_train[[split_var]]),size = 1)

                        if((sum(x_train[g_node$obs_train,split_var]==split_var_sampled_rule)>=node_min_size) & (sum(x_train[g_node$obs_train,split_var]!=split_var_sampled_rule)>=node_min_size)){
                                good_tree_index <- 1
                        } else {
                                split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                                if(length(split_var_candidates)==0){
                                        return(tree) # There are no valid candidates for this node
                                }
                        }


                }       else       {

                        # Avoiding invalid max.
                        if((length(x_train[g_node$obs_train,split_var])-node_min_size)<1){
                                return(tree)
                        }

                        # Getting the min and maximum observed value within the terminal node
                        min_node_obs <- sort(x_train[g_node$obs_train,split_var])[node_min_size]
                        max_node_obs <- sort(x_train[g_node$obs_train,split_var])[length(x_train[g_node$obs_train,split_var])-node_min_size]


                        # Getting the column from xcut
                        xcut_valid <- xcut[which(xcut[,split_var]>=min_node_obs & xcut[,split_var]<=max_node_obs),split_var]


                        # No valid tree found
                        if(length(xcut_valid) == 0 ){

                                split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                                if(length(split_var_candidates)==0){
                                        return(tree) # There are no valid candidates for this node
                                }

                        } else {
                                good_tree_index <- 1

                                # Getting only unique values of xcut_valid
                                xcut_valid <- unique(xcut_valid)
                                split_var_sampled_rule <- sample(xcut_valid,size = 1)
                        }
                }# End for no categorical variables
        }

        # Creating the left and the right nodes
        max_index <- max(get_indexes(tree))
        max_tree_size <- length(tree)

        if(split_var %in% cat_var){
                # Creating the vector of the splitting rules for the categorical
                left_train_id <- g_node$obs_train[which(x_train[g_node$obs_train,split_var]==split_var_sampled_rule)]
                right_train_id <- g_node$obs_train[which(x_train[g_node$obs_train,split_var]!=split_var_sampled_rule)]

                left_test_id <- g_node$obs_test[which(x_test[g_node$obs_test,split_var]==split_var_sampled_rule)]
                right_test_id <- g_node$obs_test[which(x_test[g_node$obs_test,split_var]!=split_var_sampled_rule)]

        } else {
                # Creating the vector of new train and test index
                left_train_id <- g_node$obs_train[which(x_train[g_node$obs_train,split_var]<=split_var_sampled_rule)]
                right_train_id <- g_node$obs_train[which(x_train[g_node$obs_train,split_var]>split_var_sampled_rule)]

                left_test_id <- g_node$obs_test[which(x_test[g_node$obs_test,split_var]<=split_var_sampled_rule)]
                right_test_id <- g_node$obs_test[which(x_test[g_node$obs_test,split_var]>split_var_sampled_rule)]
        }
        # No valid tree
        if((length(left_train_id) < node_min_size) || (length(right_train_id)<node_min_size)){
                return(tree)
        }

        # Creating the left node
        left_node <- list(index = max_index+1,
                          obs_train = left_train_id,
                          obs_test  = left_test_id,
                          left = NA,
                          right = NA,
                          parent = g_node$index,
                          terminal = 1,
                          nog = 0,
                          depth = g_node$depth+1,
                          var = split_var_sampled_rule,
                          var_split_rule = split_var_sampled_rule,
                          mu = 0)

        right_node <- list(index = max_index+2,
                          obs_train = right_train_id,
                          obs_test  = right_test_id,
                          left = NA,
                          right = NA,
                          parent = g_node$index,
                          terminal = 1,
                          nog = 0,
                          depth = g_node$depth+1,
                          var = split_var,
                          var_split_rule = split_var_sampled_rule,
                          mu = 0)

        # Modifying the new g_node
        new_g_node <-g_node
        new_g_node$left <- max_index+1
        new_g_node$right <- max_index+2
        new_g_node$terminal <- 0
        new_g_node$nog = 1 # It always be a a nog since is given origin to two children


        # Get nog counter ( FOR THE NEW TREE )
        nog_counter <- count_nog(tree = tree[-g_node_position]) + 1

        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood(res_vec = res_vec,node = left_node,tau = tau,tau_mu = tau_mu) + node_loglikelihood(res_vec = res_vec,node = right_node, tau = tau,tau_mu = tau_mu) - node_loglikelihood(res_vec = res_vec,node = g_node,tau = tau, tau_mu = tau_mu)

        # Calculate the transition
        transition_loglike <- log(0.3/nog_counter)-log(0.3/length(terminal_nodes)) # prob of getting from the new tree to the old (PRUNE), minus getting to the old to the new (GROW)

        # Calculate the tree prior contribution
        tree_prior <- 2*log(1-alpha*(1+(g_node$depth+1))^(-beta)) + ((-beta)*log((1+g_node$depth)) + log(alpha)) - log(1-alpha*(1+g_node$depth)^(-beta))

        log_acceptance <- tree_loglikeli+transition_loglike+tree_prior

        # Accepting the tree ornot
        if(stats::runif(n = 1)<=exp(log_acceptance)){
                # Maybe use append to make everything easier
                tree[[g_node_name]] <- new_g_node

                # Transform the parent into a nog
                if(!is.na(g_node$parent)){
                        tree[[paste0("node_",g_node$parent)]]$nog <- 0
                }


                new_nodes <- list(left_node,right_node)
                names(new_nodes) <- c(paste0("node_",c(new_nodes[[1]]$index,new_nodes[[2]]$index)))
                tree <- append(tree,new_nodes,after = g_node_position_orig)
        }


        tree_validator(tree = tree)
        return(tree)

}


# Creating a function to grow a tree
grow_rotation <- function(res_vec,
                 tree,
                 x_train,
                 x_test,
                 xcut,
                 tau,
                 tau_mu,
                 alpha,
                 beta,
                 node_min_size,
                 # Passing rotation variables
                 rotation_variables){

        # Getting the terminal nodes
        terminal_nodes <- get_terminals(tree)

        # Sampling one terminal node
        g_node_position <- sample(1:length(terminal_nodes),size = 1)
        g_node <- terminal_nodes[[g_node_position]]
        g_node_name <- names(terminal_nodes[g_node_position])
        g_node_position_orig <- which(names(tree) ==g_node_name)

        # Initializing the sample
        split_var_candidates <- 1:length(rotation_variables)
        good_tree_index <- 0


        while(good_tree_index==0){
                # Selecting a valid split
                split_var_pair <- sample(rotation_variables,size = 2)

                split_var <- sample(split_var_pair,size = 1)

                # Selecting an angle to rotate my coordinates
                theta <- stats::runif(n = 1,min = 0,max = pi)

                # Creating the rotated coordinates
                rotated_x <- tcrossprod(A(theta), x_train[,split_var_pair])
                rownames(rotated_x) <- split_var_pair


                # Getting the rotation for the test
                rotated_x_test <- tcrossprod(A(theta), x_test[,split_var_pair])
                rownames(rotated_x_test) <- split_var_pair

                if((length(rotated_x[split_var,g_node$obs_train])-node_min_size)<1){
                        return(tree)
                }
                # Getting the min and maximum observed value within the terminal node
                min_node_obs <- sort(rotated_x[split_var,g_node$obs_train])[node_min_size]
                max_node_obs <- sort(rotated_x[split_var,g_node$obs_train])[length(rotated_x[split_var,g_node$obs_train])-node_min_size]

                # Getting the x_cut matrix rotated
                xcut_rotated <- tcrossprod(A(theta), xcut[,split_var_pair])
                rownames(xcut_rotated) <- split_var_pair

                # Getting the column from xcut
                xcut_valid <- xcut_rotated[split_var,which(xcut_rotated[split_var,]>=min_node_obs & xcut_rotated[split_var,]<=max_node_obs)]


                # No valid tree found
                if(length(xcut_valid) == 0 ){

                        rotation_variables <-  rotation_variables[!(rotation_variables %in% split_var_pair)]

                        if(length(rotation_variables)==0 || length(rotation_variables)==1){
                                return(tree) # There are no valid candidates for this node
                        }

                } else {
                        good_tree_index <- 1
                }
        }
        # Sampling a x_cut_rule
        xcut_valid <- unique(xcut_valid)
        split_var_sampled_rule_rotation <- sample(xcut_valid,size = 1)

        # Creating the left and the right nodes
        max_index <- max(get_indexes(tree))
        max_tree_size <- length(tree)

        # Creating the vector of new train and test index
        left_train_id <- g_node$obs_train[which(rotated_x[split_var,g_node$obs_train]<=split_var_sampled_rule_rotation)]
        right_train_id <- g_node$obs_train[which(rotated_x[split_var,g_node$obs_train]>split_var_sampled_rule_rotation)]

        left_test_id <- g_node$obs_test[which(rotated_x_test[split_var,g_node$obs_test]<=split_var_sampled_rule_rotation)]
        right_test_id <- g_node$obs_test[which(rotated_x_test[split_var,g_node$obs_test]>split_var_sampled_rule_rotation)]

        # No valid tree
        if((length(left_train_id) < node_min_size) || (length(right_train_id)<node_min_size)){
                return(tree)
        }

        # Creating the left node
        left_node <- list(index = max_index+1,
                          obs_train = left_train_id,
                          obs_test  = left_test_id,
                          left = NA,
                          right = NA,
                          parent = g_node$index,
                          terminal = 1,
                          nog = 0,
                          depth = g_node$depth+1,
                          var = list(split_var_pair = split_var_pair, split_var = split_var , theta = theta),
                          var_split_rule = split_var_sampled_rule_rotation,
                          mu = 0)

        right_node <- list(index = max_index+2,
                           obs_train = right_train_id,
                           obs_test  = right_test_id,
                           left = NA,
                           right = NA,
                           parent = g_node$index,
                           terminal = 1,
                           nog = 0,
                           depth = g_node$depth+1,
                           var = list(split_var_pair = split_var_pair, split_var = split_var, theta = theta),
                           var_split_rule = split_var_sampled_rule_rotation,
                           mu = 0)

        # Modifying the new g_node
        new_g_node <-g_node
        new_g_node$left <- max_index+1
        new_g_node$right <- max_index+2
        new_g_node$terminal <- 0
        new_g_node$nog = 1 # It always be a a nog since is given origin to two children


        # Get nog counter ( FOR THE NEW TREE )
        nog_counter <- count_nog(tree = tree[-g_node_position]) + 1

        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood(res_vec = res_vec,node = left_node,tau = tau,tau_mu = tau_mu) + node_loglikelihood(res_vec = res_vec,node = right_node, tau = tau,tau_mu = tau_mu) - node_loglikelihood(res_vec = res_vec,node = g_node,tau = tau, tau_mu = tau_mu)

        # Calculate the transition
        transition_loglike <- log(0.3/nog_counter)-log(0.3/length(terminal_nodes)) # prob of getting from the new tree to the old (PRUNE), minus getting to the old to the new (GROW)

        # Calculate the tree prior contribution
        tree_prior <- 2*log(1-alpha*(1+(g_node$depth+1))^(-beta)) + ((-beta)*log((1+g_node$depth)) + log(alpha)) - log(1-alpha*(1+g_node$depth)^(-beta))

        log_acceptance <- tree_loglikeli+transition_loglike+tree_prior

        # Accepting the tree ornot
        if(stats::runif(n = 1)<=exp(log_acceptance)){
                # Maybe use append to make everything easier
                tree[[g_node_name]] <- new_g_node

                # Transform the parent into a nog
                if(!is.na(g_node$parent)){
                        tree[[paste0("node_",g_node$parent)]]$nog <- 0
                }


                new_nodes <- list(left_node,right_node)
                names(new_nodes) <- c(paste0("node_",c(new_nodes[[1]]$index,new_nodes[[2]]$index)))
                tree <- append(tree,new_nodes,after = g_node_position_orig)
        }

        return(tree)

}

# Pruning a tree
prune <- function(tree,
                  res_vec,
                  tau,
                  tau_mu,
                  alpha,
                  beta){

        # Getting the node
        nog_nodes <- get_nog(tree = tree)
        t_nodes <- (get_terminals(tree = tree))

        n_terminal_nodes <- length(t_nodes)
        n_nogs <- length(nog_nodes)

        # Returning the  a simple tree
        if(length(nog_nodes)==0){
                return(tree)
        }

        # Sample a node to be pruned
        nog_nodes_index <- sample(1:length(nog_nodes),size = 1)

        if(nog_nodes_index!=0 & length(nog_nodes)!=0){
                p_node <- nog_nodes[[nog_nodes_index]]
        }


        # Name node to be pruned
        p_node_name <- paste0("node_",p_node$index)

        # Getting the name of non terminals
        names_non_terminals <- names(tree[!(names(tree) %in% names(t_nodes))])
        names_non_terminals <- names_non_terminals[names_non_terminals!=p_node_name] # Removing the current pruned node
        parents_non_t_nodes <- sapply(tree[names_non_terminals],function(x){x$parent})
        left_node_name <- paste0("node_",p_node$left)
        left_node <- tree[[paste0("node_",p_node$left)]]
        right_node_name <- paste0("node_",p_node$right)
        right_node <- tree[[paste0("node_",p_node$right)]]

        # Calculating the loglikelihood of the new tree
        tree_loglikeli <- node_loglikelihood(res_vec = res_vec,node = p_node,tau = tau, tau_mu = tau_mu) -
                          node_loglikelihood(res_vec = res_vec,node = left_node,tau = tau,tau_mu = tau_mu) -
                          node_loglikelihood(res_vec = res_vec,node = right_node, tau = tau,tau_mu = tau_mu)

        # Calculate the transition
        transition_loglike <- log(0.3/(n_terminal_nodes-1))-log(0.3/length(nog_nodes)) # prob of getting from the new tree to the old (PRUNE), minus getting to the old to the new (GROW)

        # Calculate the tree prior contribution
        tree_prior <- log(1-alpha*(1+p_node$depth)^(-beta)) - ((-beta)*log((1+p_node$depth)) + log(alpha)) - 2*log(1-alpha*(1+(p_node$depth+1))^(-beta))

        log_acceptance <- tree_loglikeli+transition_loglike+tree_prior

        # Deciding weather accept or not
        if(stats::runif(n = 1) < exp(log_acceptance)){
                tree[[p_node_name]]$terminal <- 1
                tree[[p_node_name]]$left <- NA
                tree[[p_node_name]]$right <- NA
                tree[[p_node_name]]$nog <- 0

                # Removing the pruned nodes
                tree[[left_node_name]] <- NULL
                tree[[right_node_name]] <- NULL

                # Checking if AFTER pruning this node its parent become a NOG
                if(p_node_name!="node_0"){
                        new_t_nodes_names <- c(names(t_nodes),p_node_name)
                        pruned_node_parent <- tree[[paste0("node_",tree[[p_node_name]]$parent)]]
                        if((paste0("node_",pruned_node_parent$left) %in% new_t_nodes_names) & (paste0("node_",pruned_node_parent$right) %in% new_t_nodes_names)){
                                tree[[paste0("node_",tree[[p_node_name]]$parent)]]$nog <- 1
                        }
                }
        }

        return(tree)
}

# Changing the tree
change <- function(res_vec,
                   tree,
                   x_train,
                   x_test,
                   xcut,
                   tau,
                   tau_mu,
                   alpha,
                   beta,
                   node_min_size,
                   cat_var){


        # Getting the node
        nog_nodes <- get_nog(tree = tree)
        n_terminal_nodes <- length(get_terminals(tree = tree))


        # cat("NoG: ", length(nog_nodes), "   ||  Tree size", length(tree), "\n" )

        # Sample a node to be pruned
        nog_nodes_index <- sample(1:length(nog_nodes),size = 1)

        if(nog_nodes_index!=0 & length(nog_nodes)!=0){
                c_node <- nog_nodes[[nog_nodes_index]]
        }

        # Case of just one split
        if( (length(nog_nodes)==0) & (length(tree)==3)){
                c_node <- tree[[1]] ## Getting the root node
        }

        if(length(tree)<3){
                return(tree)
        }

        good_tree_index <- 0

        # Getting the name of the changed node
        split_var_candidates <- colnames(x_train)


        if(is.null(cat_var)){
                cat_var <- names(x_train)[!apply(x_train,2,function(x){length(unique(x))>round(nrow(x_train)*0.4)})]
        }

        while(good_tree_index==0){

                # Selecting a valid split
                split_var <- sample(split_var_candidates,size = 1)

                if(split_var %in% cat_var){

                        # Selecting the split rule for categorical variables
                        split_var_sampled_rule <- sample(unique(x_train[[split_var]]),size = 1)

                        if((sum(x_train[c_node$obs_train,split_var]==split_var_sampled_rule)>=node_min_size) & (sum(x_train[c_node$obs_train,split_var]!=split_var_sampled_rule)>=node_min_size)){
                                good_tree_index <- 1
                        } else {
                                split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                                if(length(split_var_candidates)==0){
                                        return(tree) # There are no valid candidates for this node
                                }
                        }

                } else {


                        # Case of invalid max
                        if((length(x_train[c_node$obs_train,split_var])-node_min_size)<1){
                                return(tree)
                        }


                        # Getting the min and maximum observed value within the terminal node
                        min_node_obs <- sort(x_train[c_node$obs_train,split_var])[node_min_size]
                        max_node_obs <- sort(x_train[c_node$obs_train,split_var])[length(x_train[c_node$obs_train,split_var])-node_min_size]

                        # Getting the column from xcut
                        xcut_valid <- xcut[which(xcut[,split_var]>=min_node_obs & xcut[,split_var]<=max_node_obs),split_var]


                        # No valid tree found
                        if(length(xcut_valid) == 0 ){

                                split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                                if(length(split_var_candidates)==0){
                                        return(tree) # There are no valid candidates for this node
                                }

                        } else {
                                # Sampling a x_cut_rule
                                # xcut_valid <- unique(xcut_valid)
                                split_var_sampled_rule <- sample(xcut_valid,size = 1)
                                good_tree_index <- 1
                        }

                }



        }


        if(split_var %in% cat_var){

                # Creating the vector of new train and test index
                left_train_id <- c_node$obs_train[which(x_train[c_node$obs_train,split_var]==split_var_sampled_rule)]
                right_train_id <- c_node$obs_train[which(x_train[c_node$obs_train,split_var]!=split_var_sampled_rule)]

                left_test_id <- c_node$obs_test[which(x_test[c_node$obs_test,split_var]==split_var_sampled_rule)]
                right_test_id <- c_node$obs_test[which(x_test[c_node$obs_test,split_var]!=split_var_sampled_rule)]

        } else {
                # Creating the vector of new train and test index
                left_train_id <- c_node$obs_train[which(x_train[c_node$obs_train,split_var]<=split_var_sampled_rule)]
                right_train_id <- c_node$obs_train[which(x_train[c_node$obs_train,split_var]>split_var_sampled_rule)]

                left_test_id <- c_node$obs_test[which(x_test[c_node$obs_test,split_var]<=split_var_sampled_rule)]
                right_test_id <- c_node$obs_test[which(x_test[c_node$obs_test,split_var]>split_var_sampled_rule)]
        }


        # Getting the left and the right
        new_left_name <- paste0("node_",c_node$left)
        new_right_name <- paste0("node_",c_node$right)

        # No valid tree
        if((length(left_train_id) < node_min_size) || (length(right_train_id)<node_min_size)){
                return(tree)
        }

        # Creating a new left node and changing it
        old_left_node <- tree[[new_left_name]]
        new_left_node <- tree[[new_left_name]]
        new_left_node$obs_train <- left_train_id
        new_left_node$obs_test <- left_test_id
        new_left_node$var <- split_var
        new_left_node$var_split_rule <- split_var_sampled_rule
        new_left_node$parent <- c_node$index

        # Creating a new right node and changing it
        old_right_node <- tree[[new_right_name]]
        new_right_node <- tree[[new_right_name]]
        new_right_node$obs_train <- right_train_id
        new_right_node$obs_test <- right_test_id
        new_right_node$var <- split_var
        new_right_node$var_split_rule <- split_var_sampled_rule
        new_right_node$parent <- c_node$index

        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood(res_vec = res_vec,node = new_left_node,tau = tau,tau_mu = tau_mu) +
                node_loglikelihood(res_vec = res_vec,node = new_right_node,tau = tau,tau_mu = tau_mu) -
                node_loglikelihood(res_vec = res_vec,node = old_left_node, tau = tau,tau_mu = tau_mu) -
                node_loglikelihood(res_vec = res_vec,node = old_right_node,tau = tau, tau_mu = tau_mu)

        log_acceptance <- tree_loglikeli




        # Accepting the tree ornot
        if(stats::runif(n = 1)<=exp(log_acceptance)){

                # Maybe use append to make everything easier
                tree[[new_left_name]] <- new_left_node
                tree[[new_right_name]] <- new_right_node

        }

        tree_validator(tree = tree)

        return(tree)


}


# Changing the tree
change_rotation <- function(res_vec,
                   tree,
                   x_train,
                   x_test,
                   xcut,
                   tau,
                   tau_mu,
                   alpha,
                   beta,
                   node_min_size,
                   rotation_variables){


        # Getting the node
        nog_nodes <- get_nog(tree = tree)
        n_terminal_nodes <- length(get_terminals(tree = tree))

        # print(" ROTATION TESTS!!!")

        # If there is only the stump
        if(length(tree)==1){
                return(tree)
        }

        # Sample a node to be pruned
        nog_nodes_index <- sample(1:length(nog_nodes),size = 1)

        if(nog_nodes_index!=0 & length(nog_nodes)!=0){
                c_node <- nog_nodes[[nog_nodes_index]]
        }


        # Case of just one split
        if( (length(nog_nodes)==0) & (length(tree)==3)){
                c_node <- tree[[1]] ## Getting the root node
        }

        if(length(tree)<3){
                return(tree)
        }

        # Case of just one split
        # if(length(nog_nodes)==0 & length(tree)==3){
        #         c_node <- tree[[1]] ## Getting the root node
        # }  else {
        #         return(tree)
        # }

        good_tree_index <- 0



        while(good_tree_index==0){

                # Selecting a valid pair split
                split_var_pair <- sample(rotation_variables,size = 2)

                # Selecting the splitted var from the par
                split_var <- sample(split_var_pair, size = 1)

                # Selecting an angle to rotate my coordinates
                theta <- stats::runif(n = 1,min = 0,max = pi)

                # Creating the rotated coordinates
                rotated_x <- tcrossprod(A(theta), x_train[,split_var_pair])
                rownames(rotated_x) <- split_var_pair

                # Getting the rotation for the test
                rotated_x_test <- tcrossprod(A(theta), x_test[,split_var_pair])
                rownames(rotated_x_test) <- split_var_pair

                # Case of invalid max
                if((length(rotated_x[split_var,c_node$obs_train])-node_min_size)<1){
                        return(tree)
                }

                # Getting the min and maximum observed value within the terminal node
                min_node_obs <- sort(rotated_x[split_var,c_node$obs_train])[node_min_size]
                max_node_obs <- sort(rotated_x[split_var,c_node$obs_train])[length(rotated_x[split_var,c_node$obs_train])-node_min_size]

                # Getting the x_cut matrix rotated
                xcut_rotated <- tcrossprod(A(theta), xcut[,split_var_pair])
                rownames(xcut_rotated) <- split_var_pair

                # Getting the column from xcut
                xcut_valid <- xcut_rotated[split_var,which(xcut_rotated[split_var,]>=min_node_obs & xcut_rotated[split_var,]<=max_node_obs)]


                # No valid tree found
                if(length(xcut_valid) == 0 ){

                        rotation_variables <-  rotation_variables[!(rotation_variables %in% split_var_pair)]

                        if(length(rotation_variables)==0 || length(rotation_variables)==1){
                                return(tree) # There are no valid candidates for this node
                        }

                } else {
                        good_tree_index <- 1
                }
        }
        # Sampling a x_cut_rule
        xcut_valid <- unique(xcut_valid)
        split_var_sampled_rule_rotation <- sample(xcut_valid,size = 1)


        # Creating the vector of new train and test index
        left_train_id <- c_node$obs_train[which(rotated_x[split_var,c_node$obs_train]<=split_var_sampled_rule_rotation)]
        right_train_id <- c_node$obs_train[which(rotated_x[split_var,c_node$obs_train]>split_var_sampled_rule_rotation)]

        left_test_id <- c_node$obs_test[which(rotated_x_test[split_var,c_node$obs_test]<=split_var_sampled_rule_rotation)]
        right_test_id <- c_node$obs_test[which(rotated_x_test[split_var,c_node$obs_test]>split_var_sampled_rule_rotation)]


        # No valid tree
        if((length(left_train_id) < node_min_size) || (length(right_train_id)<node_min_size)){
                return(tree)
        }

        # Getting the left and the right
        new_left_name <- paste0("node_",c_node$left)
        new_right_name <- paste0("node_",c_node$right)

        # Creating a new left node and changing it
        old_left_node <- tree[[new_left_name]]
        new_left_node <- tree[[new_left_name]]
        new_left_node$obs_train <- left_train_id
        new_left_node$obs_test <- left_test_id
        new_left_node$var <- list(split_var_pair = split_var_pair, split_var = split_var, theta = theta )
        new_left_node$var_split_rule <- split_var_sampled_rule_rotation
        new_left_node$parent <- c_node$index

        # Creating a new right node and changing it
        old_right_node <- tree[[new_right_name]]
        new_right_node <- tree[[new_right_name]]
        new_right_node$obs_train <- right_train_id
        new_right_node$obs_test <- right_test_id
        new_right_node$var <- list(split_var_pair = split_var_pair, split_var = split_var, theta = theta )
        new_right_node$var_split_rule <- split_var_sampled_rule_rotation
        new_right_node$parent <- c_node$index


        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood(res_vec = res_vec,node = new_left_node,tau = tau,tau_mu = tau_mu) +
                node_loglikelihood(res_vec = res_vec,node = new_right_node,tau = tau,tau_mu = tau_mu) -
                node_loglikelihood(res_vec = res_vec,node = old_left_node, tau = tau,tau_mu = tau_mu) -
                node_loglikelihood(res_vec = res_vec,node = old_right_node,tau = tau, tau_mu = tau_mu)

        log_acceptance <- tree_loglikeli


        # No valid tree
        if((length(new_left_node$obs_train) < node_min_size) || (length(new_right_node$obs_train)<node_min_size)){
                return(tree)
        }

        # print(exp(log_acceptance))
        # Accepting the tree ornot
        if(stats::runif(n = 1)<=exp(log_acceptance)){

                # Maybe use append to make everything easier
                tree[[new_left_name]] <- new_left_node
                tree[[new_right_name]] <- new_right_node

        }

        return(tree)


}


# Normalize BART function (Same way as theOdds code)
normalize_bart <- function(y) {

        # Defining the a and b
        a <- min(y)
        b <- max(y)

        # This will normalize y between -0.5 and 0.5
        y  <- (y - a)/(b - a) - 0.5
        return(y)
}

# Now a function to return everything back to the normal scale

unnormalize_bart <- function(z, a, b) {
        # Just getting back to the regular BART
        y <- (b - a) * (z + 0.5) + a
        return(y)
}

# Naive sigma_estimation
naive_sigma <- function(x,y){

        # Getting the valus from n and p
        n <- length(y)

        # Getting the value from p
        p <- ifelse(is.null(ncol(x)), 1, ncol(x))

        # Adjusting the df
        df <- data.frame(x,y)
        colnames(df)<- c(colnames(x),"y")

        # Naive lm_mod
        lm_mod <- stats::lm(formula = y ~ ., data =  df)

        # Getting sigma
        sigma <- summary(lm_mod)$sigma
        return(sigma)
}

# Update mu over trees
update_mu <- function(tree,
                      partial_residuals,
                      tau,tau_mu){

        # Getting the terminal nodes
        terminal_nodes <- get_terminals(tree = tree)

        for(i in 1:length(terminal_nodes)){

                terminal_node_name <- paste0("node_",terminal_nodes[[i]]$index)

                curr_partial_residuals <- partial_residuals[terminal_nodes[[i]]$obs_train]

                sum_r <- sum(curr_partial_residuals)
                curr_n_train <- length(curr_partial_residuals)

                tree[[terminal_node_name]]$mu = stats::rnorm(n = 1,mean = (sum_r*tau)/(tau*curr_n_train+tau_mu),sd = (curr_n_train*tau+tau_mu)^(-1/2))
        }

        return(tree)
}


# Getting prediction
getPrediction <- function(tree,x_train,x_test){

        # Getting the train and the test prediction
        train_pred <- rep(NA,nrow(x_train))
        test_pred <- rep(NA,nrow(x_test))

        # Getting the terminal nodes
        terminal_nodes <- get_terminals(tree = tree)

        if(length(terminal_nodes)<1){
                stop("Error NO terminal nodes in this tree")
        }
        for(i in 1:length(terminal_nodes)){

                train_pred[terminal_nodes[[i]]$obs_train] <- terminal_nodes[[i]]$mu
                test_pred[terminal_nodes[[i]]$obs_test] <- terminal_nodes[[i]]$mu

        }

        if(any(is.na(train_pred)) || any(is.na(test_pred))){
                stop(" Error in the TREE NA in the sets")
        }

        return(list(train_pred = train_pred,
                    test_pred = test_pred))
}

# Get the values
update_tau <- function(y,
                       y_hat,
                       a_tau,
                       d_tau){
        # Getting the values
        n <- length(y)

        return(stats::rgamma(n = 1,shape = (0.5*n+a_tau),rate = 0.5*crossprod( (y-y_hat))+d_tau ))
}

# Replicating the BART model
bart <- function(x_train,
                 y_train,
                 x_test,
                 n_tree,
                 n_mcmc,
                 n_burn,
                 n_min_size,
                 tau,
                 alpha, beta,
                 df, sigquant,
                 numcut,
                 scale_boolean = TRUE,
                 K_bart = 2){

        # Saving a_min and b_max
        a_min <- NULL
        b_max <- NULL

        a_min <- min(y_train)
        b_max <- max(y_train)

        # Cut matrix
        xcut <- matrix(NA,ncol = ncol(x_train),nrow = numcut)

        # Getting possible x values
        for(j in 1:ncol(x_train)){
                xs <- stats::quantile(x_train[ , j], type=7,
                                      probs=(0:(numcut+1))/(numcut+1))[-c(1, numcut+2)]

                xcut[,j] <-xs
        }# Error of the matrix
        if(is.null(colnames(x_train)) || is.null(colnames(x_test)) ) {
                stop("Insert a valid NAMED matrix")
        }

        if(!is.vector(y_train)) {
                stop("Insert a valid y vector")
        }

        # Scale values
        if(scale_boolean) {
                # Normalizing y
                y_scale <- normalize_bart(y = y_train)

                # Calculating \tau_{\mu} based on the scale of y
                tau_mu <- (4 * n_tree * (K_bart^2))

                # Getting the naive sigma
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the shape
                a_tau <- df/2

                # Calculating lambda
                qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
                lambda <- (nsigma*nsigma*qchi)/df
                d_tau <- (lambda*df)/2

        } else {

                # Not scaling the y
                y_scale <- y_train

                # Calculating \tau_{\mu} based on the scale of y
                # Need to change this value in case of non-scaling
                tau_mu <- (4 * n_tree * (K_bart^2))/((b_max-a_min)^2)
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the naive sigma
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the shape
                a_tau <- df/2

                # Calculating lambda
                qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
                lambda <- (nsigma*nsigma*qchi)/df
                d_tau <- (lambda*df)/2

        }

        # Defining other quantities
        n_train <- nrow(x_train)
        n_test <- nrow(x_test)
        n_post <- n_mcmc-n_burn

        # Getting the y_hat for train and test
        y_train_hat_post <- matrix(0, ncol = n_train,nrow = n_post)
        y_test_hat_post <- matrix(0, ncol = n_test,nrow = n_post)
        curr <- 0
        y_train_hat_trees <- matrix(0, nrow = n_tree, ncol = n_train)
        y_test_hat_trees <- matrix(0, nrow = n_tree, ncol = n_test)


        tau_post <- numeric()

        # Getting initial trees
        current_trees <- list()

        for(i in 1:n_tree){
                current_trees[[i]] <- new_tree(x_train = x_train,x_test = x_test)
        }


        # Setting the progress bar
        pb <- utils::txtProgressBar(min = 0,max = n_mcmc,style = 3)
        for(i in 1:n_mcmc){

                # Adding the tick
                Sys.sleep(0.1)
                utils::setTxtProgressBar(pb,i)

                for(t in 1:n_tree){

                        partial_residuals <- y_scale - colSums(y_train_hat_trees[-t,,drop = FALSE])

                        # Selecting one verb
                        verb <- sample(c("grow","prune","change"), prob = c(0.3,0.3,0.4),size = 1)

                        if(length(current_trees[[t]])==1){
                                verb <- "grow"
                        }

                        # Selecting a new tree
                        current_trees[[t]]  <- if(verb=="grow"){
                                grow(res_vec = partial_residuals,tree = current_trees[[t]],
                                     x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                     tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = 1)
                        } else if(verb=="prune"){
                                prune(tree = current_trees[[t]],res_vec = partial_residuals,
                                      tau = tau,tau_mu = tau_mu,alpha = alpha,beta = beta)
                        } else if(verb=="change"){
                                change(res_vec = partial_residuals,tree = current_trees[[t]],
                                       x_train = x_train,x_test = x_test,xcut = xcut,
                                       tau = tau,tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = 1)
                        }

                        # Updating mu
                        current_trees[[t]] <- update_mu(tree = current_trees[[t]],partial_residuals = partial_residuals,
                                                        tau = tau,tau_mu = tau_mu)

                        # Prediction aux
                        pred_obj <- getPrediction(tree = current_trees[[t]])

                        y_train_hat_trees[t,] <- pred_obj$train_pred
                        y_test_hat_trees[t,] <- pred_obj$test_pred

                }

                # Storing tau and getting new tau
                tau <- update_tau(y = y_scale,y_hat = colSums(y_train_hat_trees),a_tau = a_tau,d_tau = d_tau)

                tau_post[i] <- tau

                # Storing the posterior elements
                if(i>n_burn){
                        curr = curr + 1
                        y_train_hat_post[curr,] <- colSums(y_train_hat_trees)
                        y_test_hat_post[curr,] <- colSums(y_test_hat_post)
                }


        }


        # Adjusting tau and y_hat for the scale factor
        if(scale_boolean){
                tau_ost <- tau_post/((b_max-a_min)^2)
                y_train_hat_post <- unnormalize_bart(z = y_train_hat_post,a = a_min,b = b_max)
                y_test_hat_post <- unnormalize_bart(z = y_test_hat_post,a = a_min,b = b_max)
        }

        # Returning the posterior objets
        return(list(tau_post = tau_post,
                    y_hat_post = y_train_hat_post,
                    y_test_hat_post = y_test_hat_post,
                    last_trees = current_trees))

}


# BART validator -- get if the tree yielded valid tree
tree_validator <- function(tree){
        t_nodes <- get_terminals(tree = tree)
        all_t_index <- sort(unlist(lapply(t_nodes, function(x){x$obs_train})))

        if(length(all_t_index)!=length(tree[["node_0"]]$obs_train)){
                stop("Tree not valid")
        }

}
