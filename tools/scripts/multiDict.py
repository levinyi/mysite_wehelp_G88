def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].setdefault(key_b, []).append(val)
    else:
        thedict.update({key_a: {key_b: [val]}})
    return 

def threedimdict(thedict, key_a, key_b, key_c, value):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        if key_b in thedict[key_a]:
            thedict[key_a][key_b].update({key_c: value})
        else:
            thedict[key_a].update({key_b: {key_c: value}})
    else:
        thedict.update({key_a:{key_b: {key_c: value}}})
    return thedict

