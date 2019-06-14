def getDeltaT(intensity='all'):
    dictionary = {15 : -2.48084192117676
,					20 : -7.6866446317581
,					25 : -13.4333756989217
,					30 : -20.0729501617316
,					35 : -28.3905674923444
,					40 : -37.4514482130471
,					45 : -47.0349688770657
,					50 : -56.1191510037752
,					55 : -67.0731359805213
,					60 : -78.8296360661043
,					60 : -91.5736228900787
}
    if intensity == "all":
        return dictionary
    else:
        return dictionary[intensity]
def getConstants(intensity='all'):
    dictionary = {15 : 1.70250205341e-06,
					20 : 3.75210165491e-05,
					25 : 7.28923756371e-05,
					30 : 0.000108250702381,
					35 : 0.000143602638145,
					40 : 0.000178949568808,
					45 : 0.000214291959608,
					50 : 0.000249630006983,
					55 : 0.000284963805289,
					60 : 0.000320293403128,
					60 : 0.000355618826099}

    if intensity == "all":
        return dictionary
    else:
        return dictionary[intensity]
