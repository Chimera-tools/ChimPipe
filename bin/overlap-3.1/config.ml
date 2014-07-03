(* config.ml *)
(* module to manage the parameters input to the program, and to display help *)
(* main change wrt ~/bin/OcamlSrc/OverlapTwoSetsOfFeat code : inclusion type option added *)


let string_from_incltype incltype =
  match incltype with
      1 -> "i1"
    | 2 -> "i2"
    | _ -> "ov";;
    

(*********************************************************)
(* Structure for the parameters of the overlap program   *)
(*********************************************************)
type 'a context_t =
	{ 
	  mutable file1: string;		(* input file1 *)
 	  mutable file2: string;		(* input file2 *) 
	  mutable maxnbseq: int;                (* maximum number of genomic sequences in the input files *)
	  mutable labfeat2: string;		(* label for file2 features *) 
	  mutable mode: string;                 (* mode of the overlap output *)
	  mutable nonred: bool;                 (* whether non redundant file2 features overlapping a file1 feature will be displayed *)
	  mutable strmode: int;                 (* whether and how strand is taken into account *)
	  mutable outfile: string;              (* output file *)
	  mutable verbose: bool;                (* if we want the output to be displayed on stdout *)
	  mutable incltype: int;                (* inclusion type *)
	  mutable ucsc: bool;                   (* ucsc option, prints values of the 9th field with two double quotes 
						   and semicolon to output in ucsc *)
	  mutable inter: bool;                  (* intersection mode, in case the negative mode is used, reports intersection segments 
						   rather than file2 feature segments *)
	  mutable sorted: bool;                 (* boolean for whether file1 is already sorted according to chr, start end 
						   (sort -k1,1 -k4,4n -k5,5n) *)
	} 



(********************************************************)
(* overlap context, these are the default parameters    *)
(********************************************************)
let context = 
  {	
    file1 = "";
    file2= "";
    maxnbseq=50;
    labfeat2="feat2";
    mode="0";
    nonred=false;
    strmode=0; 
    outfile="";
    verbose=false;
    incltype=0;
    ucsc=false;
    inter=false;
    sorted=false;
  };;


let usage = 
" 
                     ********************************************          
                     *   overlap - version v3.1 (March 2010) *
                     *      Sarah Djebali, Sylvain Foissac      *
                     ********************************************

Usage: "^(Filename.basename (Sys.argv.(0)))^" file1 file2 [options] 

For each file1 feature, reports some information about the file2 features overlapping it.

** file1 and file2 must be provided in gff format.
** [options] can be:
   -f flag:       labelling flag for file2 features in the result file.
                  -> default is \"feat2\".

   -m mode:       mode defines the kind of overlap information one wants to retrieve about file2 features overlapping file1 features.
                   * mode=0: boolean overlap: 1 or 0 whether a given file1 feature is overlapped by a file2 feature.
                   * mode=1: number of file2 features overlapping a given file1 feature. 
                             Note: file2 features with identical {chr,start,stop} from N lines will always be counted N times. 
                   * mode=n where n<0: list and number of file2 features coordinates (in the form of chr_start_stop_strand) 
                             corresponding to file2 features overlapping a given file1 feature. 
                             Note: identical {chr,start,stop} from N lines will be provided N times, unless the -nr option is used.
                   * mode=n where n>=10 and n is even: list and number of values from the nth field of file2 corresponding 
                             to file2 features overlapping a given file1 feature.
                             Note: identical nth field values from N lines will be provided N times, unless the -nr option is used.
                   * mode=n1,n2,...,np where for each i, ni>=10 and ni is even: for each i, provides the list and the number of values 
                             from the nith field of file2 corresponding to file2 features overlapping a given file1 feature.
                             Note: identical nith field values from N lines will be provided N times, unless the -nr option is used.
                   -> default is 0 (faster).

   -nr:           outputs non redundant information about file2 features overlapping a given file1 feature. 
                  Is useful in combination with the three last -m modes (see -m option for more details). 
                  -> default is unset.

   -inter:        when negative mode is used, outputs intersection segments rather than file2 feature segments. 
                      
   -o outfile:    outfile is the name of the gff file the user wants the output to be provided in. 
                  Note: This file is file1 with some additional fields representing the overlap information one wants to retrieve.
                  -> default is file1_over[flag].gff.

   -v:            provides the output result in the standard output rather than in an output file.
                  -> default is unset.

   -st strmode:   defines whether and how strand is taken into account when computing overlap between two features:
                   * strmode=0: strand is not taken into account, only positions are
                   * strmode>0: only features on the same strand can overlap
                   * strmode<0: only features on different strands can overlap
                  Note: in case a non null strmode is used, unstranded features will not be overlapped by anything.
                   -> default is 0.
   
   -i incltype:   this option enables the users to retrieve inclusion information rather than general overlap information
                  (note that inclusion is a particular example of overlap).
                  There are two types of inclusion:
                   * incltype=2: file2 features included in file1 features.
                   * incltype=1: file2 features including file1 features.
                   -> default is 0 (general overlap).

   -ucsc:         format the output file in a way that complies with the ucsc browser 
		  (in order to directly load the file in the ucsc browser)
		  (namely adds two double quotes and one semi-colon to each value of a (key,value) pair).
                  -> default is unset.

   -s nbseq:      nbseq is an upper bound for the number of sequences you have in your input gff files.
                  -> default is 50.

   -so:           overlap does not require any sorting of any file, however this option enables to skip
                  the file1 sorting in case this file is already sorted according to chromosome, start, end. 
                  Note: this sorting could be performed outside overlap using the unix sort command: 
                  sort -k1,1 -k4,4n -k5,5n file1 > sortedfile1

** Please report any bug to sarahqd@gmail.com        
"

(*
  (unstranded features can match anything)
  * strmode=-1: only features on different strands can overlap (unstranded features can match anything)
*)



(***********************************************************************)
(* Read the arguments from the command line and updates the context    *)
(***********************************************************************)
let read_commandline () =
  let u = try 
      begin
	context.file1 <- Sys.argv.(1);
	context.file2 <- Sys.argv.(2);
      end
    with
      | Invalid_argument s -> Common.print_error usage
  in
  
  (* we start reading the arguments from the 3rd one since the two first are compulsory and are the input files *)
  let argnum = ref 2 and ok = ref true in

  (* This function returns the next argument. mustbe says whether this argument has to exist.
     In general mustbe is true since the arguments go by pairs *)
  let getarg mustbe =
    incr argnum; 
    if !argnum < (Array.length Sys.argv) then 
      Sys.argv.(!argnum)
    else if mustbe then 
      raise Not_found 
    else
      ""
  in
    (* Actually reading each of the arguments *)
    try 
      while (!ok) do
	match (getarg false) with
	  | "-f"  -> context.labfeat2 <- getarg true
	  | "-m"  -> context.mode <- getarg true
	  | "-nr" -> context.nonred <- true
          | "-st" -> context.strmode <- int_of_string (getarg true)  
	  | "-o"  -> context.outfile <- getarg true
	  | "-v"  -> context.verbose <- true
	  | "-s"  -> context.maxnbseq <- int_of_string (getarg true)
	  | "-i"  -> context.incltype <- int_of_string (getarg true)
	  | "-ucsc" -> context.ucsc <- true
	  | "-inter" -> context.inter <- true
	  | "-so" -> context.sorted <- true
	  | "-h"
	  | "--help" -> Printf.printf "%s\n" usage; exit 1; 
	  | ""	-> ok := false
	  | s	-> failwith s
      done;
      Common.print_log "# Command line has been read\n";
    with
      | Not_found -> Common.print_error ("Missing parameter at the end of the command line\n"^usage^"\n");
      | Failure "int_of_string" -> Common.print_error ("Syntax error (incorrect integer value)\n"^usage^"\n");
      | Failure s -> Common.print_error ("Syntax error ("^s^")\n"^usage^"\n");;


