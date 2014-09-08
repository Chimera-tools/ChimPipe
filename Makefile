SHELL := /bin/bash

#######################
#        Rules        #
#######################
all: write run clean
		
write: ChimPipe.sh
	@echo
	@echo "1)Write a temporary bash script"  
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "#!bin/bash" > temp.sh
	@echo "root=\$$(cd \$$(dirname \$$0); pwd -P )" >> temp.sh
	@echo "sed -i \"s|root=.*|root=\$$root|\" ChimPipe.sh" >> temp.sh
	@echo "done"
	
	
run:
	@echo
	@echo "2) Execute the script to add code to ChimPipe.sh to define an environmental" 
	@echo "variable with the directory where the pipeline is placed " 
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"	
	@$(shell bash temp.sh)
	@echo "done"
	
# Clean up
clean: temp.sh
	@echo
	@echo "3) Remove temporary bash script"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@rm -f temp.sh 
	@echo "done"


