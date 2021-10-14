needleman_wunsch <- function(seq1, seq2, match, mismatch, gap){
  buffer <- c(0:length(seq1)*gap, gap);
  buffer_length <- length(buffer);
  print(buffer);
  # seq1 is "above table" -> length(seq1) = number of columns
  # length(seq2) = number of rows
  r <- 1;
  c <- 1;
  # ptr <- buffer_length; # determines where next calculated value will be stored
  ptr <- 1;
  diagonal <- 1;
  above <- 2;
  right <- buffer_length;
  m <- matrix(nrow = length(seq2), ncol = length(seq1));
  while(TRUE){
    # print("row: ");
    # print(r);
    # print(seq2[r]);
    # print("column: ");
    # print(c);
    # print(seq1[c]);
    
    print("---------------");
    print("diag ind.:");
    print(diagonal);
    print("diag val.:");
    print(buffer[diagonal]);
    print("above ind.:");
    print(above);
    print("above val.:");
    print(buffer[above]);
    print("right ind.:");
    print(right);
    print("right val.:");
    print(buffer[right]);
    print("----");
    print("diag + mismatch:");
    print(buffer[diagonal] + mismatch);
    print("right + gap:");
    print(buffer[right]  + gap);
    print("above + gap:");
    print(buffer[above]  + gap);
    if(seq1[c] == seq2[r]){
      buffer[ptr] <- max(c(buffer[diagonal] + match, buffer[above] + gap, buffer[right]  + gap));
      # debug info
      m[r,c] <- buffer[ptr];
    }else{
      # print("b d: ");
      # print(buffer[diagonal]);
      # print("mismatch:");
      # print(mismatch);
      # print("r:");
      # print(r);
      # print("c");
      # print(c);
      
      buffer[ptr] <- max(c(buffer[diagonal] + mismatch, buffer[above] + gap, buffer[right]  + gap));
      
      # print("max");
      # print(buffer[ptr]);
      
      m[r,c] <- buffer[ptr];
      # print(m);
    }
    diagonal <- (diagonal + 1) %% buffer_length;
    above <- (above %% buffer_length + 1);
    right <- (right %% buffer_length + 1);
    ptr <- (ptr %% buffer_length + 1);
    c <- c + 1;
    
    # adjust pointers at the end of row and initialize gap in zeroth column
    if(c > length(seq1)){
      c <- 1;
      r <- r + 1;
      
      
      if(r > length(seq2)){
        print(m)
        if(ptr == buffer_length - 1){
          return(buffer[1:ptr]);
        }else if(ptr == buffer_length){
          return(buffer[2:ptr]);
        }else{
          return(c(buffer[ptr+2:buffer_length], buffer[1:ptr]));
        }
        
      }
      
      buffer[diagonal] <- gap * r;
      rigth <- diagonal;
      diagonal <- (diagonal%% buffer_length + 1) ;
      above <- (above %% buffer_length + 1);
      ptr <- (ptr %% buffer_length + 1);
      
    }
  }
}
seq1 = unlist(strsplit("TATGC", ""));
seq2 = unlist(strsplit("AGTA", ""));

# seq1 = unlist(strsplit("CGT", ""));
# seq2 = unlist(strsplit("AC", ""));
needleman_wunsch(seq1, seq2, 2, -1, -2);
