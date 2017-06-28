# Main Makefile

TOP=.

include $(TOP)/buildconfig.mk

SUBDIRS=$(TOP)/src

all:
	@for i in ${SUBDIRS}; do ( \
    	 cd $$i; \
    	 $(MAKE)) ; done

clean:
	@for i in ${SUBDIRS}; do ( \
	 cd $$i; \
	 $(MAKE) clean ) ; done
