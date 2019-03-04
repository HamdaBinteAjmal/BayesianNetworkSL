# 2007 Rainer Opgen-Rhein
# License: GNU GPL 2 or later


makedotfile_VAR <- function(results, filename)
{
   f <- file(filename, "w")  # open an output file connection

   #require(graph)
 

   # thresholds for line width and coloring
   cutoff <- quantile(abs(results[,1]), c(0.2, 0.8)) 

   cat("Digraph G {\n", file=f)
   cat("K=\"0.3\";\n", file=f)
   cat("ratio=\"0.75\";\n", file=f)
   #cat("label = \"\\n\\nNetwork\\ndrawn by FDP\";\n", file=f)
   #cat("fontsize=20;\n", file=f)
   for(i in 1:nrow(results))   
   {
#      cat(results[i,2],"->",results[i,3], "[", file=f)
       cat("\"",results[i,2],"\"->\"",results[i,3], "\"  [", file=f)

     # direction
#      if(results[i,7]=="1to2"){
       cat("dir=\"forward\",", file=f)
#      }else if (results[i,7]=="2to1"){
#        cat("dir=\"back\",", file=f)
#      }else{
#        cat("dir=\"none\",", file=f)
#      }
   
                  


      # line thickness and color depends on relative strengh	
      if (abs(results[i,1]) < cutoff[1]){    # lower 20% quantile
        if(results[i,1] < 0){
           cat("color=grey, style=\"dashed,setlinewidth(1)\"", file=f)
         }else{
           cat("color=grey, style=\"solid,setlinewidth(1)\"", file=f)
         } 
      }else if(abs(results[i,1]) < cutoff[2]){      # from 20%-80%
         if(results[i,1] < 0){
           cat("color=black, style=\"dashed,setlinewidth(1)\"", file=f)
         }else{
           cat("color=black, style=\"solid,setlinewidth(1)\"", file=f)
         }
      }else{         # top 80%-100%
         if(results[i,1] < 0){
           cat("color=black, style=\"dashed,setlinewidth(2)\"", file=f)
         }else{
           cat("color=black, style=\"solid,setlinewidth(2)\"", file=f)
         }
      }

      cat("];\n", file=f)
   }
   cat("}\n", file=f)
   close(f)
}
