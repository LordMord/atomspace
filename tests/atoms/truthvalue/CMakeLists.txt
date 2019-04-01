
LINK_LIBRARIES(
	smob
	atomspace
	value
	distvalue
)

ADD_CXXTEST(SimpleTruthValueUTest)
ADD_CXXTEST(TVUTest)
ADD_CXXTEST(TVMergeUTest)
ADD_CXXTEST(EvidenceCountTruthValueUTest)
# ADD_CXXTEST(IndefiniteTruthValueUTest)
