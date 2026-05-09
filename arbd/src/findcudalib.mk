################################################################################
#
# Copyright 1993-2013 NVIDIA Corporation.  All rights reserved.
#
# NOTICE TO USER:   
#
# This source code is subject to NVIDIA ownership rights under U.S. and 
# international Copyright laws.  
#
# NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
# CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
# IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
# REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
# MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.   
# IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
# OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE 
# OR PERFORMANCE OF THIS SOURCE CODE.  
#
# U.S. Government End Users.  This source code is a "commercial item" as 
# that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of 
# "commercial computer software" and "commercial computer software 
# documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995) 
# and is provided to the U.S. Government only as a commercial end item.  
# Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
# 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
# source code with only those rights set forth herein.
#
################################################################################
#
#  findcudalib.mk is used to find the locations for CUDA libraries and other
#                 Unix Platforms.  This is supported Mac OS X and Linux.
#
################################################################################

## Find Location of most recent CUDA Toolkit
ifeq (,$(CUDA_PATH))
    CUDA_PATH := $(shell echo $(wildcard /usr/local/cuda*) | tr ' ' '\n' | tail -n1)
    ifeq (,$(CUDA_PATH))
        $(error Could not CUDA_PATH. Please pass as follows: $(MAKE) CUDA_PATH=/path/to/cuda)
    endif
    $(info Using CUDA_PATH=$(CUDA_PATH))
endif

# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr "[:lower:]" "[:upper:]")
OSLOWER = $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/" -e "s/armv7l/32/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# Determine OS platform and unix distribution
ifeq ("$(OSLOWER)","linux")
   # first search lsb_release
   DISTRO  = $(shell lsb_release -i -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
   DISTVER = $(shell lsb_release -r -s 2>/dev/null)
   ifeq ("$(DISTRO)",'') 
     # second search and parse /etc/issue
     DISTRO = $(shell more /etc/issue | awk '{print $$1}' | sed '1!d' | sed -e "/^$$/d" 2>/dev/null | tr "[:upper:]" "[:lower:]")
     DISTVER= $(shell more /etc/issue | awk '{print $$2}' | sed '1!d' 2>/dev/null
   endif
   ifeq ("$(DISTRO)",'') 
     # third, we can search in /etc/os-release or /etc/{distro}-release
     DISTRO = $(shell awk '/ID/' /etc/*-release | sed 's/ID=//' | grep -v "VERSION" | grep -v "ID" | grep -v "DISTRIB")
     DISTVER= $(shell awk '/DISTRIB_RELEASE/' /etc/*-release | sed 's/DISTRIB_RELEASE=//' | grep -v "DISTRIB_RELEASE")
   endif
endif

# search at Darwin (unix based info)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))
ifneq ($(DARWIN),)
   SNOWLEOPARD = $(strip $(findstring 10.6, $(shell egrep "<string>10\.6" /System/Library/CoreServices/SystemVersion.plist)))
   LION        = $(strip $(findstring 10.7, $(shell egrep "<string>10\.7" /System/Library/CoreServices/SystemVersion.plist)))
   MOUNTAIN    = $(strip $(findstring 10.8, $(shell egrep "<string>10\.8" /System/Library/CoreServices/SystemVersion.plist)))
   MAVERICKS   = $(strip $(findstring 10.9, $(shell egrep "<string>10\.9" /System/Library/CoreServices/SystemVersion.plist)))
   MAVERICKS   = $(strip $(findstring 10.9, $(shell egrep "<string>10\.9" /System/Library/CoreServices/SystemVersion.plist)))
endif 

# Common binaries
GCC   ?= g++
CLANG ?= /usr/bin/clang++

ifeq ("$(OSUPPER)","LINUX")
	CC=$(GCC)
else
    # for some newer versions of XCode, CLANG is the default compiler, so we need to include this
    ifeq ($(shell expr `xcodebuild -version | grep -i xcode | awk '{print $$2}' | cut -d'.' -f1` \>= 5),1)
        CC = $(CLANG)
		CC_FLAGS += -stdlib=libstdc++
		NV_FLAGS += -Xcompiler -arch -Xcompiler x86_64 -Xcompiler -stdlib=libstdc++
    endif
endif
NVCC ?= $(CUDA_PATH)/bin/nvcc -ccbin $(CC)

# Take command line flags that override any of these settings
ifeq ($(i386),1)
	OS_SIZE = 32
	OS_ARCH = i686
endif
ifeq ($(x86_64),1)
	OS_SIZE = 64
	OS_ARCH = x86_64
endif
ifeq ($(ARMv7),1)
	OS_SIZE = 32
	OS_ARCH = armv7l
endif
