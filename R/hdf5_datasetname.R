datasetname_to_df <- function (datasetname) {
    words_list <- stringr::str_split(datasetname, '_')
    
    res_list <- lapply(words_list, function (words) {
        res <- list(diagram = words[1], quarks = words[2])
        ops_list <- stringr::str_split(words[3:length(words)], '\\.')
        ops_parsed <- lapply(ops_list, function (ops) {
            p_str <- ops[1]
            d_str <- ops[2]
            g_str <- ops[3]
            
            p_parsed <- as.numeric(stringr::str_match(p_str, 'p(-?\\d)(-?\\d)(-?\\d)')[1, 2:4])
            d_parsed <- as.numeric(stringr::str_match(d_str, 'd(-?\\d)(-?\\d)(-?\\d)')[1, 2:4])
            g_parsed <- as.numeric(stringr::str_match(g_str, 'g(-?\\d)')[1, 2])
            
            list(p = p_parsed, d = d_parsed, g = g_parsed)
        })
        
        for (j in 1:4) {
            for (d in 1:3) {
                res[[sprintf('p%i%s', j, c('x', 'y', 'z')[d])]] <- ops_parsed[[j]]$p[d]
                res[[sprintf('d%i%s', j, c('x', 'y', 'z')[d])]] <- ops_parsed[[j]]$d[d]
            }
            res[[sprintf('g%i', j)]] <- ops_parsed[[j]]$g
        }
        
        res
    })
    
    res_t <- lapply(paramvalf::list_transpose(res_list), unlist)
    df <- as.data.frame(res_t)
    df$datasetname <- datasetname
    return (df)
}

df_to_datasetname <- function (df) {
    # We need to figure out how many operators are involved. This is the number
    # of columns starting with “p”.
    matches <- stringr::str_match(colnames(df), 'p\\dx')[, 1]
    op_num <- sum(!is.na(matches))
    
    words <- sapply(1:op_num, function (op_id) {
        sprintf('p%d%d%d.d%d%d%d.g%d',
                df[[sprintf('p%dx', op_id)]],
                df[[sprintf('p%dy', op_id)]],
                df[[sprintf('p%dz', op_id)]],
                df[[sprintf('d%dx', op_id)]],
                df[[sprintf('d%dy', op_id)]],
                df[[sprintf('d%dz', op_id)]],
                df[[sprintf('g%d', op_id)]])
    })
    joined <- apply(words, 1, function (row) paste(row, collapse = '_'))
    datasetname <- sprintf('%s_%s_%s', df$diagram, df$quarks, joined)
}

p_q_to_df_2 <- function (p, q) {
    pq <- p + q
    
    data.frame(p1x = pq[1],
               p1y = pq[2],
               p1z = pq[3],
               p2x = -q[1],
               p2y = -q[2],
               p2z = -q[3])
}

p_q_to_df_4 <- function (p, q1, q2) {
    pq1 <- p + q1
    pq2 <- p + q2
    
    data.frame(p1x = pq1[1],
               p1y = pq1[2],
               p1z = pq1[3],
               p2x = -q1[1],
               p2y = -q1[2],
               p2z = -q1[3],
               p3x = pq2[1],
               p3y = pq2[2],
               p3z = pq2[3],
               p4x = -q2[1],
               p4y = -q2[2],
               p4z = -q2[3])
}