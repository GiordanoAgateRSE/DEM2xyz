#-------------------------------------------------------------------------------
# "DEM2xyz v.3.0" (DEM manager tool)
# Copyright 2016-2022 (RSE SpA)
# "DEM2xyz v.3.0" authors and email contact are provided on the documentation 
# file.
# This file is part of DEM2xyz v.3.0 .
# DEM2xyz v.3.0 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# DEM2xyz is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with DEM2xyz. If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Description: DEM2xyz Makefile.
#-------------------------------------------------------------------------------

# Variables to be modified: start

#    CODE VERSION
#       Suggested values:
#          " vx_UserInitials_CommitData"
VERSION = v3_

#    COMPILER
#       Suggested values:
#          " ifort"
#          " gfortran"
COMPILER = gfortran

#    EXECUTION
#       Suggested values:
#          " optimized"
#          " debug"
EXECUTION = debug

# Variables to be modified: end

# IDs of the above choices reported in the name of the executable code
ifeq ($(COMPILER), ifort)
   ID_COMPILER = 1
   else ifeq ($(COMPILER), gfortran)
      ID_COMPILER = 2
endif
ifeq ($(EXECUTION), optimized)
   ID_EXECUTION = 1
   else ifeq ($(EXECUTION), debug)
      ID_EXECUTION = 2
endif

# COMPILATION_FLAGS
COMPILATION_FLAGS = -O1
ifeq ($(EXECUTION), debug)
   ifeq ($(COMPILER), gfortran)
      COMPILATION_FLAGS = -g -O0 -fbacktrace -C -Wall -fbounds-check
      else ifeq ($(COMPILER), ifort)
         COMPILATION_FLAGS = -g -O0 -traceback -C -check bounds -check noarg_temp_created -debug all -fpe0
   endif
endif

EXE_DIR = ../exe/

SOURCES = \
$(wildcard ./Geometry/*.f90) \
$(wildcard ./Topography/*.f90)

MAIN_FILE_ROOT = ./Main/DEM2xyz

EXE_CODE = \
$(EXE_DIR)DEM2xyz_$(VERSION)_$(ID_COMPILER)$(ID_EXECUTION)

OBJECTS = $(patsubst %.f90,%.o,$(SOURCES))

# Targets: start
#    Correct sequence to compile+link from a Linux terminal: 
#       "make"

all: touch compile link clean 

compile: $(SOURCES)

$(SOURCES): %.f90: %.o
	$(COMPILER) $@ -o $< $(COMPILATION_FLAGS) -c

link: $(OBJECTS)
	touch $(EXE_CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 $^ -o $(EXE_CODE).x $(COMPILATION_FLAGS)
	rm -f $(MAIN_FILE_ROOT).o

clean: $(OBJECTS)
	rm -f $^
	rm -f *.mod

touch:
	touch $(OBJECTS)

remove: 
	rm -f $(EXE_CODE).x

# Targets: end
