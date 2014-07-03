(* output.ml *)


open Feature


(* print takes as input:
   - u which is the ucsc boolean saying whether we want the output to be formatted as ucsc 
   - o which is the output channel corresponding to the file where we need to write file1 with additional info 
   - f1 which is current file1 feature with additional overlap information that we want to print
   and prints in o the feature f1, in the ucsc manner if u is true, and only if the feature is not null.
   Note: ucsc manner means that we add two double quotes and a semi-colon to all values of the 9th field of the
   output that do not already contain some.
*)
let print o u f1 = 
  if(Feature.isnull f1) then
    ()
  else
    Feature.print_gff_flex o u f1;;

  


 
