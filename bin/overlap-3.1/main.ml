(* main.ml: for any feature of the first file outputs either whether or not it is overlapped by
   a feature of the second file, or by how many features of the second file, and which of them,
   it is overlapped.
   The improvement with regards to the previous overlap version (with Incl) is that a very big file1
   will not cause any ram problem since file1 is first unix sorted and then read line by line
 *)

open Common
open Config
open Feature
open Input
open Output
open Printf



(* this is the most general overlap function that given the array of file2 features corresponding
   to the sequence of the current file1 feature (line1) and a current position in this array,
   outputs the list of file2 features overlapping the current file1 feature and updates the 
   current position in this array.
   This function is meant to be used by all the other more restrictive overlap functions.
   The restriction could be due to two main things:
   - strand mode (only features on same or diff strands can overlap)
   - inclusion mode (only features included or including a file1 feature will be output.
   Then some printing restrictions can come from three things:
   - redundancy mode (the user only wants non redundant info)
   - intersection mode (the user wants the intersection segments between file1 and file2 features, only with negative mode)
   - main mode (the user wants partial information like boolean or number of file2 features overlap)
   At the end the final printing will also be different according to three things:
   - the output channel
   - the ucsc mode
   - the feature2 label
*)
let overlap_general featarr2 currpos2 line1 = 
  let nbfeat2 = Array.length featarr2 and feat2list = ref [] in
    (* while current file2 features lies upstream of current file1 feature, go forward in file2 *)
    while (((!currpos2)<nbfeat2) && ((Feature.stop (featarr2.(!currpos2)))<(Feature.start line1))) do
      incr currpos2;
    done;
    (* the fact that we want to remember the first file2 feature overlapped by the current file1 feature
       (associated to the fact that we have to go backward in file2 file afterwards) is only for the    
       computation of the overlap with the next file1 feature, not for the current one.
       (see comments in ~/bin/OcamlSrc/OverlapTwoSetsOfFeat/overlap_better_all.awk)
    *)
    let firstcurrpos2 = !currpos2 in
      (* if the current file2 feature actually overlaps the current file1 feature then we position at the first  
	 file2 features that is in this case (I am now wondering whether this first while is useful here
	 but keep it for now since it cannot harm, in the worst case we go backward one step
	 to go forward one step just afterwards...)
      *)
      if (((!currpos2)<nbfeat2)&&(Feature.start (featarr2.(!currpos2)))<=(Feature.stop line1)) then
	begin
	  while ((!currpos2)>=0&&(Common.foverlap (Feature.interval line1) (Feature.interval featarr2.(!currpos2)))) do
	    decr currpos2;
	  done;
	  incr currpos2;
	  (* and from this first file2 feature in this case we go forward in file2 in order to get
	     all possible file2 features overlapping the current file1 feature.
	     Note that the overlap condition is kept separate from the start-stop condition
	     since there could be a mixture of successive file2 features, with some that overlap
	     the current file1 feature and others not, even in between the ones that do,
	     see examples in ~/bin/OcamlSrc/Overlap_Test_on_Examples/README.sh
	  *)
	  while (((!currpos2)<nbfeat2)&&(Feature.start (featarr2.(!currpos2)))<=(Feature.stop line1)) do
	    if (Common.foverlap (Feature.interval line1) (Feature.interval featarr2.(!currpos2))) then
	      begin
		feat2list:=(featarr2.(!currpos2))::(!feat2list);
	      end;
	    incr currpos2;
	  done; (* end of while (((!currpos2)<nbfeat2)&&(Feature.start (featarr2.(!currpos2)))<=(Feature.stop line1)) *)
	  (* here we readjust the current file2 position in order to correctly make the computation of overlap with the next file1 feature *)
	  currpos2:=firstcurrpos2;
	end;  (* end of if (((!currpos2)<nbfeat2)&&(Feature.start (featarr2.(!currpos2)))<=(Feature.stop line1)) *)
      List.rev (!feat2list);;


(* The restrict_feat2_list_strand will restrict a given list of feature 2 overlapping a given feature 1
   according to the strand restriction option (passed by the user).
   Here we assume that if the wanted strand mode is not null then the file 1 feature is either of the +
   or of the - strand (filter has been done before = in overlap_main).
*)
let restrict_feat2_list_strand strand feat1 feat2list =
  let str1 = Feature.str feat1 in
    if (strand=0) then
      feat2list
    else
      begin
	if (strand>0) then
	  begin
	    List.filter (Feature.issamestr str1) feat2list
	  end
	else
	  (* here strand < 0 *)
	  begin
	    List.filter (Feature.isdiffstr str1) feat2list
	  end
      end;;


(* The restrict_feat2_list_strand will restrict a given list of feature 2 overlapping a given feature 1
   according to the inclusion restriction option (passed by the user). If the inclusion mode is neither
   1 nor 2, then we are in the general overlap mode.
*)
let restrict_feat2_list_inclusion inclusion feat1 feat2list =
  if (inclusion==2) then
    (* here we only want the feature 2 that are included in the feat1 *)
    begin
      List.filter (fun f2 -> Feature.inclusion f2 feat1) feat2list
    end
  else
    begin
      if (inclusion==1) then
	(* here we only want the features 2 that include the feat1 *)
	begin
	  List.filter (Feature.inclusion feat1) feat2list
	end
      else
	(* here general overlap mode, so we keep everything in feat2list *)
	feat2list
    end;;

  

(* The restrict_feat2_list will restrict the complete list of feature 2 overlapping a given feature 1
   according to the main restriction options which are the strand and the inclusion (passed by the user).
   Here we assume that if the wanted strand mode is not null then the file 1 feature is either of the +
   or of the - strand (filter has been done before = in overlap_main).
*)
let restrict_feat2_list strand inclusion feat1 feat2list =
  restrict_feat2_list_inclusion inclusion feat1 (restrict_feat2_list_strand strand feat1 feat2list);;

  




(****************** 5 most atomic functions *************************************)

(***** The 5 most atomic overlap functions (5 for 5 general modes: 0, 1, nofld, neg, listfld) *****)
(***** are all taking the same last 6 arguments:                                              *****)
(***** - incltype: the inclusion tye (needed for the making the new keys),                    *****)
(***** - outchan: the output channel,                                                         *****)
(***** - ucsc: the ucsc mode,                                                                 *****)
(***** - labfeat2: the feature 2 label (those 3 are for output formatting,                    *****)
(***** - feat2list: the list of file2 features overlapping the current file1 feature          *****)
(***** - feat1: the current file1 feature.                                                    *****)
(****  All are printing the current file1 line with overlap info                              *****)
(***** Then we have the three following functions that have 1 or more inner arguments:        *****)
(***** overlap_number_list nofield nr                                                         *****)
(***** overlap_number_list_feat2_coord nr inter                                               *****)
(***** overlap_numbers_lists lvalues nr                                                       *****)




(******************************* Overlap - BOOL mode *************************************)
let overlap_bool incltype outchan ucsc labfeat2 feat2list feat1 = 
  let inclstring = Config.string_from_incltype incltype in
  let string = (inclstring)^("_")^(labfeat2)^(":") in
  let booloverlap = ((List.length feat2list)!=0) in
    Output.print outchan ucsc (Feature.setattlist (List.append (Feature.attlist feat1) [(string, (Common.bool_to_string booloverlap))]) feat1);;
    


(******************************* Overlap - NUMBER mode *************************************)
let overlap_number incltype outchan ucsc labfeat2 feat2list feat1 =
  let inclstring = Config.string_from_incltype incltype in
  let string = ("nb_")^(inclstring)^("_")^(labfeat2)^(":") in
  let nbover = List.length feat2list in
    Output.print outchan ucsc (Feature.setattlist (List.append (Feature.attlist feat1) [(string, (string_of_int (nbover)))]) feat1);;



(******************************* Overlap - NUMBER_LIST mode *************************************)
let overlap_number_list nofield nr incltype outchan ucsc labfeat2 feat2list feat1 = 
  let inclstring = Config.string_from_incltype incltype in
  let first_string = ("nb_")^(inclstring)^("_")^(labfeat2)^(":") and second_string = ("list_")^(labfeat2)^(":") in
    (* 
       - f2toadd is the current file2 field to add,
       - att_list_from_feat2 is the final attribute list to add to feat1 
    *)
  let f2toadd= ref "" and att_list_from_feat2 = ref [] in
    begin
      try
	List.iter
	  (fun feat2 -> 
	    f2toadd:=(snd (List.nth (Feature.attlist feat2) ((nofield-10)/2)));
	    att_list_from_feat2:=(Common.insert_right_place nr (!f2toadd) (!att_list_from_feat2) String.compare)
	  )
	  feat2list 
      with
	  (* a failure would come from List.nth, meaning that the field no nofield does not exist *)
	| Failure s -> Common.print_error "The specified file2 field is empty"
	| _ -> Common.print_error "Problem in overlap_number_list, when making the list of attributes"
    end;
    (* seg fault around here so I put a tail recursive function instead of List.append *)
    try 
      Output.print outchan ucsc (Feature.setattlist (List.rev (((second_string,Common.list_to_string (!att_list_from_feat2)))::(first_string,(string_of_int (List.length (!att_list_from_feat2))))::(List.rev (Feature.attlist feat1)))) feat1)
    with
      | _ -> Common.print_error "Problem in overlap_number_list, when printing an output line"
    (* before we had
       Output.print outchan ucsc (Feature.setattlist (List.append (Feature.attlist feat1) [(first_string,(string_of_int (List.length (!att_list_from_feat2)))); (second_string,Common.list_to_string (!att_list_from_feat2))]) feat1);;
    *)



(******************************* Overlap - NUMBER_LIST_FEAT2_COORD mode ***********************************)
let overlap_number_list_feat2_coord nr inter incltype outchan ucsc labfeat2 feat2list feat1 = 
  let inclstring = Config.string_from_incltype incltype in
  let seq = try (Feature.seq feat1) with Invalid_argument _ -> "" in
  let first_string = ("nb_")^(inclstring)^("_")^(labfeat2)^(":") and second_string = ("list_")^(labfeat2)^(":") in
    (* 
       - begs, ends, strs are the coordinates of the current file2 feature
       - f2toadd is the current segment to add
       - att_list_from_feat2 is the final attribute list to add to feat1 
    *)
  let begs = ref "" and ends = ref "" and strs = ref "" and f2toadd= ref "" in
  let att_list_from_feat2 = ref [] in
    List.iter
      (fun feat2 ->
	if (not inter) then
	  begin
	    begs:=string_of_int (Feature.start feat2);
	    ends:=string_of_int (Feature.stop feat2);
	    strs:=strand_to_string (Feature.str feat2);
	  end
	else
	  begin
	    (* need to compute the intersection between file1 feature and current file 2 feature *)
	    begs:=string_of_int (fst (Feature.intersection feat1 feat2));
	    ends:=string_of_int (snd (Feature.intersection feat1 feat2));
	    strs:=".";
	  end;
	f2toadd:=((seq)^("_")^(!begs)^"_"^(!ends)^"_"^(!strs));
	att_list_from_feat2:=(Common.insert_right_place nr (!f2toadd) (!att_list_from_feat2) String.compare)
      )
      feat2list;
    Output.print outchan ucsc (Feature.setattlist (List.append (Feature.attlist feat1) [(first_string,(string_of_int (List.length (!att_list_from_feat2)))); (second_string,Common.list_to_string (!att_list_from_feat2))]) feat1);;



(******************************* Overlap - NUMBERS_LISTS mode *************************************)
let overlap_numbers_lists lvalues nr incltype outchan ucsc labfeat2 feat2list feat1 =  
  (* here we need several first_string and several second_string = as many as elements in lvalues *)
  let inclstring = Config.string_from_incltype incltype in
  let first_string_list = List.map (fun v -> ("nb_")^(inclstring)^("_")^(labfeat2)^("_")^(string_of_int v)^(":")) lvalues in
  let second_string_list = List.map (fun v -> ("list_")^(labfeat2)^("_")^(string_of_int v)^(":")) lvalues in
    (* we need as many f2toadd as elements in lvalues *)
  let f2toadd_list = List.map (fun v -> ref "") lvalues in
    (* we need as many att_list_from_feat2 as elements in lvalues *)
  let att_list_from_feat2_list = List.map (fun v -> ref []) lvalues in
  let reffeat1 = ref feat1 in
    begin
      try
	List.iter
	  (fun feat2 ->
	    (* here we update the list of current elements to add.
	       when we had only one field the command was
	       f2toadd:=(snd (List.nth (Feature.attlist feat2) ((nofield-10)/2)));
	    *)
	    Array.iteri 
	      (
		fun i v -> let r = (List.nth f2toadd_list i) in 
			     r:=(snd (List.nth (Feature.attlist feat2) ((v-10)/2)));
	      )
	      (Array.of_list lvalues);
	    
	    (*
	      here we update the list of list of elements.
	      when we had only one field the command was the following:
	      att_list_from_feat2:=(Common.insert_right_place nr (!f2toadd) (!att_list_from_feat2) String.compare)
	    *)
	    Array.iteri 
	      (
		fun i listref -> listref:=(Common.insert_right_place nr (!(List.nth f2toadd_list i)) (!listref) String.compare)
	      )
	      (Array.of_list att_list_from_feat2_list)
	  )
	  feat2list
      with
	| Failure s -> Common.print_error "The specified file2 field is empty"
    end;
    (* We update the attribute list of the current file1 feature.
       here we have as many attributes to add as elements in lvalues 
       when we had only one field the command was
       Output.print outchan ucsc (Feature.setattlist (List.append (Feature.attlist feat1) [(first_string,(string_of_int (List.length (!att_list_from_feat2)))); (second_string,Common.list_to_string (!att_list_from_feat2))]) feat1;
    *)
    Array.iteri
      (
	(* for each requested number of field in file2, we add a list of 2 pairs to the attribute list 
	   of the current file1 feature, representing the values of the requested fields in file2
	   of the file2 features overlapping the current file1 feature
	*)
	fun i v -> 
	  (
	     reffeat1:=(Feature.setattlist 
			   (
			     List.append 
			       (Feature.attlist (!reffeat1))
			       [
				 ((List.nth first_string_list i),(string_of_int (List.length (!(List.nth att_list_from_feat2_list i))))); 
				 ((List.nth second_string_list i),Common.list_to_string (!(List.nth att_list_from_feat2_list i)))
			       ]
			   ) (!reffeat1))))
      (Array.of_list lvalues);
    Output.print outchan ucsc (!reffeat1);;
  







(**********************************************  MAIN SUBFUNCTION OF MAIN FUNCTION  *************************************)

(*
  This function determines which is the final atomic overlap function we will use
  according to the three following main parameters provided by the user:
  - the first parameter is the overlap mode,
  - the second is the redundancy mode, 
  - the third is the intersection mode,
*)
let overlap_wanted mode redund inter =
  try 
    (
      let intmode = int_of_string mode in  (* before we were matching (intmode,redund,strand,incltype) *)
	match (intmode,redund) with
	    (0,_) -> overlap_bool 
	  | (1,_) -> overlap_number 
	  | (n,r) when (n>=10&&(n mod 2)==0) -> overlap_number_list n r
	  | (n,r) when n<0 -> overlap_number_list_feat2_coord r inter 
	  | _ -> failwith "with the -m n option and when n is an integer, n must be 0, 1, strictly negative, or greater than 10 and even\n"
    ) 
  with
      (* a failure would come from the expression 
	 let intmode = int_of_string mode 
	 and would mean that the mode is not a simple integer. In this case we have to check 
	 wether it is a lits of integers separated by commas, and that each of them is >=10 and even. 
      *)
    | Failure _ -> try 
	  (
	    let lval = get_list_of_gff_value_field mode in
	    let lvalint = List.map (int_of_string) lval in
	      overlap_numbers_lists lvalint redund
	  ) 
      with
	  (* a failure would come from get_list_of_gff_value_field and would mean that
	     one of the elements of the comma separated string is not an integer.
	  *)
	| Failure _ -> failwith "with the -m n option and when n is not an integer, n should be a list of integers greater than 10 and even, separated by commas\n";;
  


(**********************************************  MAIN FUNCTION *************************************)

(* Starts from two files and outputs one file, the first one, but to which it concatenates
   the information of whether each feature is overlapped by a feature of the second file
   or by how many features of the second file it is overlapped (and their list). 
   A main improvement eregarding the memory used by the program is not to store the whole file1
   but instead read file1 line by line.
*)
let overlap_main () =
  (* Read the arguments on the command line *)
  read_commandline ();
  
  (* Read file2 and store its info in a feature array hashtable h2, in an ORDERED way *)
  let h2 = make_feat_arr_hash context.maxnbseq context.file2 in
  let u = Common.print_log "# hashtable of feature arrays made for file2\n" in
    
  (* Makes a hashtable with same keys as h2, ie file1 sequence names, and as keys the 
       current index position in the corresponding feature array, ie at the begining
       initialized to ref 0 and then will be updated when we go along sorted file1.
  *)
  let h2indices = Hashtbl.create (Hashtbl.length h2) in
  let u = Hashtbl.iter (fun seq x -> Hashtbl.add h2indices seq (ref 0)) h2 in
  let u = Common.print_log "# I have created the hashtable of indices for file2\n" in
  
  (* Sort file1 according to chr, start, end and put the result in an intermediate file 
     tmp_sorted_file1 placed by default in the /tmp directory since a /tmp directory always 
     exists on all systems, but it may be interesting to add an option to enable the use 
     to specify his own tmp directory (in case of huge file for example). This argument
     would then need to be passed to the make_temp_file_name function.
  *)
  let tmp_sorted_file1 = Common.make_sorted_temp_file_if_user_wants context.sorted context.file1 in
  let u= Common.print_log ("# I have treated (sorted and put in temp file) file1 according to what the user wants\n") in

  (* we first open the input and output channels (Note: determine the output chanel according to 
     context.verbose and context.outfile 
  *)
  let inchan1 = open_in tmp_sorted_file1 and stop = ref false and i = ref 0 in
  let outchan = 
    if (context.verbose) then
      stdout
    else
      begin
	try
	  (open_out context.outfile) 
	with
	    (* maybe we should use _incl2 or _incl1 in case inclusion information is asked for? *)
	  | Sys_error s -> open_out ((context.file1)^("_over")^(context.labfeat2)^(".gff"))
      end
  in
  
    (* then we read the temporary sorted file1 line by line (in order to avoid the RAM issue)  
       and for each line we apply the overlap algorithm depending on what kind of overlap the user wants.
       Note: takes into account context.verbose and context.ucsc to know where and how to write.
    *)
    while (not !stop) do
      (
	try
	  let currline1 = (split '\t' (input_line inchan1)) in
	  let currgffline1 = record_of_line_gff_flex currline1 in
	  let currseq1 = Feature.seq currgffline1 and currstr1 = Feature.str currgffline1 in
	    incr i;
	    (* 
	       in case the user specifies a non null strand mode and if the current file1 feature is neither on
	       the + nor on the - strand, then the current file1 feature will not be overlapped by anything
	       so we output line1 with default added information (according to parameters)
	    *)
	    if((context.strmode!=0)&&(not ((currstr1=Forward)||(currstr1=Reverse)))) then
	      ((overlap_wanted context.mode context.nonred context.inter) 
		  context.incltype outchan context.ucsc context.labfeat2 [] currgffline1)
	    else
	      begin
		let completefeat2list = overlap_general 
		  (try Hashtbl.find h2 currseq1 with Not_found -> [||]) 
		  (try Hashtbl.find h2indices currseq1 with Not_found -> (ref 0)) 
		  currgffline1 in
		let feat2list_restricted = restrict_feat2_list context.strmode context.incltype currgffline1 completefeat2list in
		  (* overlap_wanted itself writes the current temporary sorted file1 line with the info added, 
		     in order to avoid the RAM issue. In case the sequence of the current file1 line is absent 
		     in file2, an empty table is passed for h2 
		  *)
		let u= (overlap_wanted context.mode context.nonred context.inter) 
		  context.incltype outchan context.ucsc context.labfeat2 feat2list_restricted currgffline1 in
		  ()
		  (* Common.print_log (("line ")^(string_of_int (!i))^(" of sorted file 1 treated\n")) *)
	      end
	with
	  | End_of_file -> stop:=true
	  | Failure "int_of_string" -> Common.print_pb ("# Conversion problem from string to integer when reading file1\n")
	  | Stack_overflow | Out_of_memory -> Common.print_log ("# Stack overflow or Out of memory problem when computing the overlap\n")
      )
    done;  
    let u = Common.print_log ("# Overlap did its work ! Now writing output file\n") in
    let u = Common.remove_sorted_temp_file_if_user_wants context.sorted tmp_sorted_file1 in
    let u = Common.print_log ("# I have removed the temporary sorted file if it has been created\n") in
      flush outchan;;


overlap_main ();;

