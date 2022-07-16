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
        poo = line.split("=")
        #if (line.find("zmconv_ke") !=-1):
        if (poo[0].strip() == parm):
            print(poo)
            poo[1] = "  "+ valu  +" \n"
            #poo="\\".join( poopoo )
            print(poo)
            line="=".join(poo)

        fex.write(line)
    
    fin.close()
    fex.close()

    cmd = "mv "+filo+" "+fili
    sp.run( cmd, shell=True)
