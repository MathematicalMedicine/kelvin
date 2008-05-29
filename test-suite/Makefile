# Makefile to perform kelvin tests

TESTS = \
	TP_DT_LD \
	TP_DT_LE \
	TP_DT_MM_LD \
	TP_DT_MM_LE \
	TP_QT_LE \
	MP_DT \
	MP_QT

.PHONY : all $(TESTS)

all : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile test

