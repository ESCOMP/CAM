import subprocess as sp

# which namelist parameter to change

def nmled(fili,parm,valu):
    #which file
    #fili  = "atm_in"
    filo  = fili+"_edit"


    fin   = open( fili ,"r")
    fex   = open( filo ,"w")
    linin = fin.readlines()
    for line in linin:
        spl = line.split("=")
        #if (line.find("zmconv_ke") !=-1):
        if (spl[0].strip() == parm):
            print(spl)
            spl[1] = "  "+ valu  +" \n"
            print(spl)
            line="=".join(spl)

        fex.write(line)
    
    fin.close()
    fex.close()

    cmd = "mv "+filo+" "+fili
    sp.run( cmd, shell=True)

def nmlread(fili,parm):

    valu=-99999
    fin   = open( fili ,"r")
    linin = fin.readlines()
    for line in linin:
        spl = line.split("=")
        if (spl[0].strip() == parm):
            valu=spl[1]
    
    fin.close()

    return valu
