(*******************************************************************************)
(* Ce fichier contient des petites fonctions qui sont utilisees un peu partout *)
(*******************************************************************************)


let identity x y = 
  if (x=y) then
    0
  else
    1

(*************************)
(* Printing functions    *)
(*************************)

(* Indique si on est en mode verbeux *)
let verbose = ref true;;

(* Permet d'afficher des logs seulement si le flag verbose est positionne *)
let print_log s =
	if !verbose then output_string stderr s;
	flush stderr;;

(* Permet d'afficher un message d'erreur et de quitter *)
let print_error s =
	output_string stderr s;
	flush stderr;
	exit 1;;

let print_pb s =
	output_string stderr ((s)^"\n");
	flush stderr;;

(* in order to be sure that two users do not use the program at the same second on the machine
   we could add a random number to index the temp file. Also it would be good to be able to place
   the file elsewhere than in the /tmp directory of the system, in case of huge file for example.
   In this case the user could specify a place and make_temp_file_name will take another parameter.
   In fact rather than a random number I will add the process id.
*)
let make_temp_file_name file =
  let time = Unix.localtime (Unix.time ()) and pid = Unix.getpid () in 
    (("/tmp/")^(string_of_int (1+time.Unix.tm_mon))^("_")^(string_of_int time.Unix.tm_mday)^("_")^(string_of_int time.Unix.tm_hour)^"_"^(string_of_int time.Unix.tm_min)^("_")^(string_of_int time.Unix.tm_sec)^("_")^(string_of_int pid));;

let make_sorted_temp_file_if_user_wants sorted file =
  let tmpfile=make_temp_file_name file in
    if (not sorted) then
      begin
	Unix.system (("sort -k1,1 -k4,4n -k5,5n ")^(file)^(" > "^(tmpfile)));
	tmpfile;
      end
    else
      file;;

let remove_sorted_temp_file_if_user_wants sorted file =
  if (not sorted) then
    begin
      Unix.system (("rm ")^(file));
      ();
    end;;
      


(*******************************************)
(* Functions on integer intervals          *)
(*******************************************)
(* general overlap *)
let foverlap i1 i2 =
  let (beg1,end1)=i1 and (beg2,end2)=i2 in
  (end1>=beg2)&&(beg1<=end2);;

(* general inclusion (not strict one) of the first interval in the second one *)
let finclusion i1 i2 =
  let (beg1,end1)=i1 and (beg2,end2)=i2 in
    (beg1>=beg2)&&(end1<=end2);;

(* strict inclusion of the first interval in the second one *)
let fstrictinclusion i1 i2 =
  let (beg1,end1)=i1 and (beg2,end2)=i2 in
    (beg1>beg2)&&(end1<end2);;



(*******************************************)
(* Functions on strings                    *)
(*******************************************)

let bool_to_string b =
  match b with
      true -> "1"
    | false -> "0";;

(* Extract le nom de la proteine, de la sequence ... en enlevant les commentaires
   ie: ce qu'il y a apres un espace ou un | pipe *)
let extract_name s deb =
	let t1 = try
		String.index s ' '
	with
		Not_found -> (String.length s) in
	let t2 = try
		min t1 (String.index s '|')
	with
		Not_found -> t1 in
	String.sub s deb (t2-deb);;


(* suffix s i retourne le suffixe d'une chaine s a partir de la position i *)
let suffix s i = 
	try
		String.sub s i ((String.length s)-i)
	with
		Invalid_argument "String.sub" -> "";;



(* split c s découpe la chaine de caractères s selon le caractere c *)
let split c s = 
	let rec split_from n = 
		try
			let p = String.index_from s n c
			in (String.sub s n (p-n)) :: (split_from (p+1))
		with Not_found -> [ suffix s n ] 
	in if s="" then [] else split_from 0 ;;

(* rev_split fait comme split mais renvoie la liste a l'envers -> recursivite terminale *)
(* rev_split enleve le dernier bout de la chaine s'il est egal a "" *)
let rev_split c s =
	let rec rev_split_from n acc =
		try
			let p = String.index_from s n c
			in rev_split_from (p+1) ((String.sub s n (p-n))::acc)
		with Not_found -> match suffix s n with 
							| "" -> acc
							| p -> p::acc
	in if s="" then [] else rev_split_from 0 [];;


(* clean_end_string eliminates all space characters present at the end of a string.
   This is useful to clean the last field of a gff file, since it can contain
   such spaces at the end. Note that it cannot contain any \t character since
   these characters are separating the 8 first fields only. 
*)
let clean_end_string s =
  let n = String.length s in
  let i= ref (n-1) in
  let unclean=ref(s.[!i]=' ') in
    while ((!unclean) && (!i>=1)) do
      decr i;
      unclean:=(s.[!i]=' ');
    done;
    if((!i=0)&&(!unclean)) then
      ""
    else
      String.sub s 0 (!i+1);;
    
(* clean_beg_string eliminates all space characters present at the begining of a string.
   This is useful to clean the last field of a gff file, since it can contain
   such spaces. Note that it cannot contain any \t character since
   these characters are separating the 8 first fields only. 
*)
let clean_beg_string s =
  let n = String.length s in
  let i= ref 0 in
  let unclean=ref(s.[!i]=' ') in
    while ((!unclean) && (!i<=n-2)) do
      incr i;
      unclean:=(s.[!i]=' ');
    done;
    if((!i=(n-1))&&(!unclean)) then
      ""
    else
      String.sub s (!i) (n-(!i));;
    


(****************************)
(* Functions on lists *)
(****************************)

let rec intervals l =
  match l with
    |[] -> [];
    |[_] -> [];
    |t1::t2::q -> (t1,t2)::(intervals (t2::q))


let in_interval n (i,j) =
  (n>=i) && (n<=j)

(* Cette fontion renvoie le dernier element d'une liste *)
let rec last = function
	| [] -> failwith "Liste vide"
	| [x] -> x
	| _::l -> last l;;

(* on veut tri_fusionner deux listes, selon un ordre sur leurs elements donné par comp *)
(*let rec tri_fusion comp l1 l2 =
  match l1,l2 with
    |(l,[]) -> l;
    |([],l) -> l;
    |(t1::q1,t2::q2) -> if(comp t1 t2 <= 0) then 
	t1::(tri_fusion comp q1 (t2::q2)) 
      else 
	t2::(tri_fusion comp (t1::q1) q2);;*)
let tri_fusion = List.merge;;



(* ordonne prend une liste de listes d'hsps et renvoie une liste d'hsps contenant toutes les hsps de toutes 
   les sous-listes de la liste donnee en paramètre, mais triées selon leur debut dans le génameique croissant *)
(*let rec flat_and_order comp = function
	| [] -> [];
	| [l1] -> l1;
	| l1::l2::ql -> flat_and_order comp ((tri_fusion comp l1 l2)::ql);;*)
let flat_and_order comp l = List.fold_left (fun a b -> tri_fusion comp a b) [] l;;


(* Enlève la redondance selon comp dans une liste supposée triée selon comp *)
let rec remove_redund comp = function
	| [] -> []
	| [t] -> [t]
	| t1::t2::q -> if (comp t1 t2) = 0 then remove_redund comp (t2::q) else t1::(remove_redund comp (t2::q));;


(* insert_right_place is a function that places an object o in a list lo ordered
   according to the function comp in the right place according to comp.
   nr is a boolean which says whether we want to insert o without or with redundancy
   in the list lo. This function is not tail recursive and may thus cause
   stack overflow problems.
*)
let rec insert_right_place nr o lo comp = 
  match lo with
    | [] -> [o]
    | t::q -> 
	if (comp o t) < 0 then 
	  o::lo
	else
	  begin
	    if (comp o t) > 0 then 
	      t::(insert_right_place nr o q comp)
	    else 
	      begin
		if nr then  (* here we ask to insert o without redundancy and o is equal to t so we do not insert it *)
		  t::q
		else  (* here we insert it *)
		  t::o::q
	      end
	  end ;;


(* retire la première occurence d'un élément dans une liste *)
let rec removelt h = function
	| [] -> []
	| t::q -> if t = h then q else t::(removelt h q);;

 
 
(* idem mais avec une liste d'elements à pb à retirer de l *)
(*let rec remove_all lpb l = 
	match lpb with
		| [] -> l
		| tpb::qpb -> if ((removelt tpb l)!=l) then 
						remove_all qpb (removelt tpb l)
					else
						match l with
						| [] -> []
						| t::q -> t::(remove_all qpb q);;*)

(* Supprime tous les elements a probleme d'une liste *)
let rec remove_all lpb l = 
	match lpb with
		| [] -> l
		| tpb::qpb -> remove_all qpb (removelt tpb l);;


let rec listmissnb ltrie = function
	| [] -> []
	| t::q -> if not (List.mem t ltrie) then t::(listmissnb ltrie q) else listmissnb ltrie q;;



let rec spread elt = function
	| [] -> []
	| t::q -> (elt::t)::(spread elt q);;


(* takes as input a list of strings and builds a string from them
   by putting "," between the different elements.
   Note: at the end of the output string there will be a comma
*)
let rec list_to_string = function
  | [] -> "";
  | t::q -> (t)^(",")^(list_to_string q);;


(******************************)
(* Fonctions sur les tableaux *)
(******************************)


(* trouve_index permet de trouver l'indice (ref : 1) d'un objet dans un tableau *)
let trouve_index tp p =
	let i = ref 0 and trouve = ref false and nbp = Array.length tp in
	while (!i<nbp && (not !trouve)) do
		trouve := (p = tp.(!i));
		incr i;
	done;
	!i-1;;




let inv_comp a b = -(Pervasives.compare a b);;


