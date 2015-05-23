SHELL := /bin/bash

#######################
#        Rules        #
#######################
all: write run clean

write:		
	@echo
	@echo "Setting up ChimPipe V0.9.0 ..."
	@echo "#!bin/bash" > temp.sh
	@echo "root=\$$(cd \$$(dirname \$$0); pwd -P )" >> temp.sh
	@echo "sed -i \"s|rootDir=.*|rootDir=\$$root|\" ChimPipe.sh" >> temp.sh
	@echo "sed -i \"s|rootDir=.*|rootDir=\$$root|\" src/bash/*.sh" >> temp.sh	
run:
	@$(shell bash temp.sh)

# Clean up
clean: temp.sh
	@rm -f temp.sh 
	@echo "done!"
	@echo

