# Makefile to perform kelvin tests

TESTS = \
	SA_DT \
	SA_QT \
	SA_QT_LC2 \
	TP_DT_LD \
	TP_DT_LE \
	TP_DT_MM_LD \
	TP_DT_MM_LE \
	TP_QT_LE

.PHONY : test $(TESTS)

test : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile test

