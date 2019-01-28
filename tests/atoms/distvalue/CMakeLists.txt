
LINK_LIBRARIES(
	atomspace
	value
	distvalue
	atomutils
)

ADD_CXXTEST(CHistTest)
ADD_CXXTEST(IntervalTest)
ADD_CXXTEST(HistogramTest)
ADD_CXXTEST(DVTest)
ADD_CXXTEST(CDVTest)
ADD_CXXTEST(DVFormulaTest)
