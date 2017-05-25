__author__ = 'mjohnpayne'

import string

text = "g fmnc wms bgblr rpylqjyrc gr zw fylb. rfyrq ufyr amknsrcpq ypc dmp. bmgle gr gl zw fylb gq glcddgagclr ylb rfyr'q ufw rfgq rcvr gq qm jmle. sqgle qrpgle.kyicrpylq() gq pcamkkclbcb. lmu ynnjw ml rfc spj."

alphabet = string.ascii_lowercase

def sub(i):
    newstring = ''
    for j in i:
        if j in alphabet[:-2]:
            j = alphabet[alphabet.index(j)+2]
            newstring += j
        elif j =='y':
            newstring += 'a'
        elif j == 'z':
            newstring += 'b'
        else:
            newstring += j
    return newstring

print sub(text)

print sub("http://www.pythonchallenge.com/pc/def/map.html")

