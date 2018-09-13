
import matplotlib.pyplot as plt
import numpy as np

####################### functions ###############
def centerPlot(cellCenterX, cellCenterY) : 

    plt.scatter(cellCenterX, cellCenterY, color="k")

    plt.xlabel('X (micro meter)')
    plt.ylabel('Y (micro meter)')
    plt.title('locations of cell cetner')
    plt.grid(True)
    plt.savefig("test.png")
    plt.show() 


####################### input parameter data ###############
nPouchCells=29  
nPeripCells=6
nBcCells=5 # each side

lumen=2 


wPouch=2 
hPouchA=10
hPouchB=10 

bNodePouch=14  #number of basal node
aNodePouch=14 #number of apical nodes
lNodePouchA=aNodePouch*5
lNodePouchB=aNodePouch*5

wPerip=10 +0.92783  
hPerip=2 

bNodePerip=5*14  #number of basal node
aNodePerip=5*14 #number of apical nodes
lNodePerip=14

cellTotalNode=168

domainMinX=-1.394 
cellsGapX=0.329  
domainMinY=14.3  
domainMaxY=domainMinY+hPouchA+lumen+hPerip  

bNodeBc=12  #number of basal node
aNodeBc=12 #number of apical nodes
lNodeBc=12

buffer_ECM=0.525
numECMBCNodes=100 ## for one side 

rBc=1  ## radius of boundary cells initialized with circular shape

################## calculation of coordinate of cell centers ######################
cellCenterPouchX=[domainMinX+cellsGapX+wPouch/2] 
cellCenterPouchY=[domainMinY+hPouchA/2] 
for x in range(nPouchCells-1):
    cellCenterPouchX.append ( cellCenterPouchX[x]+wPouch + cellsGapX)
    cellCenterPouchY.append ( cellCenterPouchY[x]) 


cellCenterPeripX =[domainMinX+cellsGapX+wPerip/2]  
cellCenterPeripY=[domainMinY+hPouchA+lumen+hPerip/2] 
for x in range (nPeripCells-1) :
    cellCenterPeripX.append ( cellCenterPeripX[x]+wPerip + cellsGapX)   
    cellCenterPeripY.append ( cellCenterPeripY[x] )  


radiusBc=(cellCenterPeripY[0]-cellCenterPouchY[0])/2
centerRightBcX=cellCenterPouchX[nPouchCells-1] + wPouch/2 + cellsGapX
centerRightBcY=0.5* ( cellCenterPeripY[0]+cellCenterPouchY[0] ) 


cellCenterBcX=[0 for x in range (nBcCells)]
cellCenterBcY=[0 for x in range (nBcCells)]

correctFactorA=1.2
correctFactorB=np.pi*0.1

for x in range (nBcCells) : 
    cellCenterBcX[x]= centerRightBcX+ radiusBc*np.sin ( (x+1)*np.pi*correctFactorA/(nBcCells+1) -correctFactorB)
    cellCenterBcY[x]= centerRightBcY- radiusBc*np.cos ( (x+1)*np.pi*correctFactorA/(nBcCells+1) -correctFactorB)


centerLeftBcX=cellCenterPouchX[0] - wPouch/2-cellsGapX
centerLeftBcY=centerRightBcY


for x in range (nBcCells) : 
    cellCenterBcX.append ( centerLeftBcX- radiusBc*np.sin ((x+1)*np.pi*correctFactorA/(nBcCells+1)-correctFactorB) )
    cellCenterBcY.append ( centerLeftBcY- radiusBc*np.cos ((x+1)*np.pi*correctFactorA/(nBcCells+1)-correctFactorB) )
    

######################### unify and sort cell centers locations ############
cellCenterX=[]
cellCenterY=[]
cellType=[]

for x in range (nPouchCells) : 
    cellCenterX.append (cellCenterPouchX[x] )  
    cellCenterY.append (cellCenterPouchY[x] )
    cellType.append ('pouch')  


for x in range (nBcCells) : 
    cellCenterX.append ( cellCenterBcX[x] )  
    cellCenterY.append ( cellCenterBcY[x] )
    cellType.append ('bc')


for x in range (1,nPeripCells+1) :
    cellCenterX.append ( cellCenterPeripX[nPeripCells-x] )  
    cellCenterY.append ( cellCenterPeripY[nPeripCells-x] ) 
    cellType.append ('peri') 

for x in range (nBcCells,2*nBcCells) : 
    cellCenterX.append ( cellCenterBcX[x] )  
    cellCenterY.append ( cellCenterBcY[x] )
    cellType.append ('bc')


####################### plot cell centers locations #####################
centerPlot(cellCenterX, cellCenterY) 


####################### Write cell centers locations as an output file #####################

fileM = open('coordinate_Cell16.txt','w') 

numCells= nPouchCells + nPeripCells + 2*nBcCells 
fileM.write('{}\n'.format(numCells))
for i in range ( numCells ) :
    #fileM.write(str(i) +"," +str (xC[i][j])   )  #,yC[i][j],typeC[i][j]) 
    fileM.write('{0:.3f}'.format(cellCenterX[i]))
    fileM.write(' ')
    fileM.write('{0:.3f}'.format(cellCenterY[i]))
    fileM.write(' ')

    fileM.write('{0:.3f}'.format(0.000))
    fileM.write(' ')

    fileM.write('{}\n'.format(cellType[i]))


fileM.close()



########################### start finding the membrane nodes ##################################################


xC=[]  ## for k =0 to create the list
yC=[]  ## for k =0 to create the list
typeC=[]

############ finding membrane nodes location of pouch cells #######################
for k in range ( nPouchCells ) :
    lastpoint=0 
    if k==0 :
        lNodePouchB=14
        lNodePouchA=42
        hPouchB=2
        hPouchA=6
        dh=2
        
    elif  k==(nPouchCells-1) :
        lNodePouchB=42
        lNodePouchA=14
        hPouchB=6
        hPouchA=2
        dh=2
    elif k==1 :
        lNodePouchB=42
        lNodePouchA=70
        hPouchB=6
        hPouchA=10
        dh=2
    elif k==(nPouchCells-2) :
        lNodePouchB=70
        lNodePouchA=42
        hPouchB=10
        hPouchA=6
        dh=2
    else:
        lNodePouchB=70
        lNodePouchA=70
        hPouchB=10
        hPouchA=10
        dh=0

    cellTotalNodePouch=lNodePouchB+lNodePouchA+bNodePouch+aNodePouch
    xC.append ( [0 for x in range ( cellTotalNodePouch)] )
    yC.append ( [0 for x in range ( cellTotalNodePouch)] )
    typeC.append ( ['notAssigned1' for x in range ( cellTotalNodePouch)] ) 

    # from center to top equal to H/2
    for i  in range (int(lNodePouchA/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA'
        xC[k][j]=cellCenterX[k]+wPouch/2 
        yC[k][j]=cellCenterY[k]+i/(lNodePouchA/2)*hPouchA/2  


    lastpoint=lastpoint+int (lNodePouchA/2)  
    for i  in range (aNodePouch) :
        j=lastpoint+i 
        typeC[k][j]='apical1'
        xC[k][j]=cellCenterX[k]+wPouch/2- i/aNodePouch*wPouch 
        if k==0 or k==1 :
            yC[k][j]=cellCenterY[k]+hPouchA/2 - dh*(i)/aNodePouch 
        elif k==(nPouchCells-2) or k==(nPouchCells-1) :
            yC[k][j]=cellCenterY[k]+hPouchA/2 + dh*i/aNodePouch  
        else:
            yC[k][j]=cellCenterY[k]+hPouchA/2 


    lastpoint=lastpoint+ aNodePouch  

    for i  in range (lNodePouchB) :
        j=lastpoint+i 
        typeC[k][j]='lateralB'
        xC[k][j]=cellCenterX[k]-wPouch/2 
        yC[k][j]=cellCenterY[k]+hPouchB/2  - i/lNodePouchB*hPouchB 


    lastpoint=lastpoint+ lNodePouchB 
    for i  in range (bNodePouch) :
        j=lastpoint+i 
        typeC[k][j]='basal1'
        xC[k][j]=cellCenterX[k]-wPouch/2+ i/bNodePouch*wPouch

        if k==0 or k==1 :
            yC[k][j]=cellCenterY[k]-hPouchB/2 -dh*i/bNodePouch
        elif k==(nPouchCells-2) or k==(nPouchCells-1) :
            yC[k][j]=cellCenterY[k]-hPouchB/2 +dh*(i)/bNodePouch
        else:
            yC[k][j]=cellCenterY[k]-hPouchB/2 


    lastpoint=lastpoint+ bNodePouch 
    for i  in range (int(lNodePouchA/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA'
        xC[k][j]=cellCenterX[k]+wPouch/2 
        yC[k][j]=cellCenterY[k]-hPouchA/2+ i/(lNodePouchA/2)*hPouchA/2  

############ finding membrane nodes location of right hand side BC cells #######################

for k in range ( nPouchCells,nPouchCells+int(nBcCells) ) :
    lastpoint=0

    cellTotalNodeBc=2*lNodeBc+bNodeBc+aNodeBc
    typeC.append ( ['notAssigned1' for x in range ( cellTotalNodeBc)] ) 
    xC.append ( [0 for x in range ( cellTotalNodeBc)] )
    yC.append ( [0 for x in range ( cellTotalNodeBc)] ) 

    # from center to top equal to H/2
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA'
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/(lNodeBc/2)*np.pi*0.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/(lNodeBc/2)*np.pi*0.25)


    lastpoint=lastpoint+int (lNodeBc/2)  
    for i  in range (int(aNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='apical1'
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/aNodeBc*np.pi*0.5+np.pi*0.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/aNodeBc*np.pi*0.5+np.pi*0.25)

    lastpoint=lastpoint+int (aNodeBc)  

    for i  in range (int(lNodeBc))  :
        j=lastpoint+i
        typeC[k][j]='lateralB'
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/lNodeBc*np.pi*0.5+np.pi*0.75)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/lNodeBc*np.pi*0.5+np.pi*0.75)


    lastpoint=lastpoint+int (lNodeBc)  
    for i  in range (int(bNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='basal1'
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/bNodeBc*np.pi*0.5+np.pi*1.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/bNodeBc*np.pi*0.5+np.pi*1.25)


    lastpoint=lastpoint+int (bNodeBc)  
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA'
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/(lNodeBc/2)*np.pi*0.25+np.pi*1.75)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/(lNodeBc/2)*np.pi*0.25+np.pi*1.75)


############ finding membrane nodes location of peripodial cells #######################

for k in range ( nPouchCells+int(nBcCells),nPouchCells+int(nBcCells)+nPeripCells ) :

    lastpoint=0
    cellTotalNodePerip=2*lNodePerip+bNodePerip+aNodePerip
    typeC.append ( ['notAssigned1' for x in range ( cellTotalNodePerip)] )
    xC.append ( [0 for x in range ( cellTotalNodePerip)] )
    yC.append ( [0 for x in range ( cellTotalNodePerip)] ) 
    # from center to top equal to H/2
    for i  in range (int(lNodePerip/2))  :
        j=lastpoint+i
        typeC[k][j]='lateralA'
        xC[k][j]=cellCenterX[k]+wPerip/2 
        yC[k][j]=cellCenterY[k]+i/(lNodePerip/2)*hPerip/2  


    lastpoint=lastpoint+int (lNodePerip/2)  
    for i  in range (bNodePerip) :
        j=lastpoint+i 
        typeC[k][j]='basal1'
        xC[k][j]=cellCenterX[k]+wPerip/2- i/bNodePerip*wPerip
        yC[k][j]=cellCenterY[k]+hPerip/2  


    lastpoint=lastpoint+ bNodePerip  
    for i  in range (lNodePerip) :
        j=lastpoint+i
        typeC[k][j]='lateralB' 
        xC[k][j]=cellCenterX[k]-wPerip/2 
        yC[k][j]=cellCenterY[k]+hPerip/2  - i/lNodePerip*hPerip 
    

    lastpoint=lastpoint+ lNodePerip 
    for i  in range (aNodePerip) :
        j=lastpoint+i 
        typeC[k][j]='apical1' 
        xC[k][j]=cellCenterX[k]-wPerip/2+ i/aNodePerip*wPerip
        yC[k][j]=cellCenterY[k]-hPerip/2  


    lastpoint=lastpoint+ aNodePerip 
    for i  in range (int(lNodePerip/2))  :
        j=lastpoint+i
        typeC[k][j]='lateralA' 
        xC[k][j]=cellCenterX[k]+wPerip/2 
        yC[k][j]=cellCenterY[k]-hPerip/2+ i/(lNodePerip/2)*hPerip/2  

############ finding location of membrane nodes of left hand side BC cells #######################

for k in range ( nPouchCells+nPeripCells+int(nBcCells),nPouchCells+nPeripCells+int(2*nBcCells) ) :

    lastpoint=0
    cellTotalNodeBc=2*lNodeBc+bNodeBc+aNodeBc
    typeC.append ( ['notAssigned1' for x in range ( cellTotalNodeBc)] )
    xC.append ( [0 for x in range ( cellTotalNodeBc)] )
    yC.append ( [0 for x in range ( cellTotalNodeBc)] ) 
    # from center to top equal to H/2
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA' 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/(lNodeBc/2)*np.pi*0.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/(lNodeBc/2)*np.pi*0.25)


    lastpoint=lastpoint+int (lNodeBc/2)  
    for i  in range (int(aNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='apical1'
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/aNodeBc*np.pi*0.5+np.pi*0.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/aNodeBc*np.pi*0.5+np.pi*0.25)


    lastpoint=lastpoint+int (aNodeBc)  
    for i  in range (int(lNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='lateralB' 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/lNodeBc*np.pi*0.5+np.pi*0.75)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/lNodeBc*np.pi*0.5+np.pi*0.75)


    lastpoint=lastpoint+int (lNodeBc)  
    for i  in range (int(bNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='basal1'
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/bNodeBc*np.pi*0.5+np.pi*1.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/bNodeBc*np.pi*0.5+np.pi*1.25)


    lastpoint=lastpoint+int (bNodeBc)  
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA' 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/(lNodeBc/2)*np.pi*0.25+np.pi*1.75)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/(lNodeBc/2)*np.pi*0.25+np.pi*1.75)


######################## plot membrane nodes locations #########################
for k in range ( numCells ) :   
    plt.scatter(xC[k], yC[k])
    plt.gca().set_aspect('equal', adjustable='box')
#plt.xlabel('X (micro meter)')
#plt.ylabel('Y (micro meter)')
#plt.title('locations of cell cetner')
#plt.grid(True)
#plt.savefig("test.png")
#plt.show() 

#centerPlot(xC[0], yC[0]) 


######################## write membrane nodes location and type #########################
numNode=[0 for x in range (numCells)]
#num_Node=[ 0 for x in range (numCells) ]
totalNode=0 
for k in range ( numCells ) :
    numNode [k]= len (xC[k]) 
    totalNode=totalNode+numNode [k]

fileM = open('coordinate_Membrane3.txt','w') 
 
#fileM.write('cellID,x coordinate, y coordinate, nodeType\n') 
for i in range ( numCells ) :
    for j in range (numNode [i]):
        #fileM.write(str(i) +"," +str (xC[i][j])   )  #,yC[i][j],typeC[i][j]) 
        fileM.write('{}'.format(i))
        fileM.write('')
        fileM.write('{0:.4f}'.format(xC[i][j]))
        fileM.write('')
        fileM.write('{0:.4f}'.format(yC[i][j]))
        fileM.write(' ')
        fileM.write('{}\n'.format(typeC[i][j]))

fileM.close() 


########################### start finding the ECM nodes location and type ##################################################
xECM=[]
yECM=[]
eCMType=[]

################# calculation of coordinates of pouch side ECM ####################
i=0
for k in range ( numCells ) :
    for j in range (numNode [k]):
        if typeC[k][j]=='basal1':
            if  cellType[k] =='pouch': 
                xECM.append(xC[k][j]) 
                yECM.append(yC[k][j]-buffer_ECM) 
                eCMType.append ('excm')
                i=i+1

################# calculation of coordinates of ECM nodes which are neighbor with right hand side BC cells ####################
enlargeR=rBc+buffer_ECM
for k in range ( numECMBCNodes) :
    xECM.append(centerRightBcX +(radiusBc+enlargeR)*np.sin ( (k+1)*np.pi/(numECMBCNodes+1) ))
    yECM.append(centerRightBcY- (radiusBc+enlargeR)*np.cos ( (k+1)*np.pi/(numECMBCNodes+1) ))
    eCMType.append ('bc2')
    i=i+1

################# calculation of coordinates of peripodial side ECM ####################
for k in range ( numCells ) :
    for j in range (numNode [k]):
         if typeC[k][j]=='basal1':
            if  cellType[k] =='peri': 
                xECM.append(xC[k][j])
                yECM.append(yC[k][j]+buffer_ECM) 
                eCMType.append ('perip')
                i=i+1
            
################# calculation of coordinates of ECM nodes which are neighbor with left hand side BC cells ####################
for k in range ( numECMBCNodes) :
    xECM.append(  centerLeftBcX- (radiusBc+enlargeR)*np.sin ((k+1)*np.pi/(numECMBCNodes+1)) )
    yECM.append(  centerLeftBcY- (radiusBc+enlargeR)*np.cos ((k+1)*np.pi/(numECMBCNodes+1)) ) 
    eCMType.append ('bc2')
    i=i+1

numECMnodes=i

########################### plotting ECM nodes ##################################################
plt.scatter(xECM, yECM, color="k")
plt.show()                


########################### writing as an output file ECM nodes locations and type ##################################################
fileM = open('coordinate_ECM16.txt','w') 
 
#fileM.write('cellID,x coordinate, y coordinate, nodeType\n')
fileM.write('{}\n'.format(numECMnodes)) 
for k in range ( numECMnodes ) :
        #fileM.write(str(i) +"," +str (xC[i][j])   )  #,yC[i][j],typeC[i][j]) 
        fileM.write('{0:.4f}'.format(xECM[k]))
        fileM.write(' ')
        fileM.write('{0:.4f}'.format(yECM[k]))
        fileM.write(' ')
        fileM.write('{}\n'.format(eCMType[k]))


fileM.close() 











