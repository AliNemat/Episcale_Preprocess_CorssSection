
import matplotlib.pyplot as plt
import numpy as np


def centerPlot(cellCenterX, cellCenterY) : 

    plt.scatter(cellCenterX, cellCenterY, color="k")

    plt.xlabel('X (micro meter)')
    plt.ylabel('Y (micro meter)')
    plt.title('locations of cell cetner')
    plt.grid(True)
    plt.savefig("test.png")
    plt.show() 



nPouchCells=29  
nPeripCells=6
nBcCells=5 # each side

lumen=2 


wPouch=2 
hPouchR=10
hPouchL=10 

bNodePouch=14  #number of basal node
aNodePouch=14 #number of apical nodes
lNodePouchR=aNodePouch*5
lNodePouchL=aNodePouch*5

wPerip=10 +0.92783  
hPerip=2 

bNodePerip=5*14  #number of basal node
aNodePerip=5*14 #number of apical nodes
lNodePerip=14

cellTotalNode=168

domainMinX=-1.394 
cellsGapX=0.329  
domainMinY=14.3  
domainMaxY=domainMinY+hPouchR+lumen+hPerip  

bNodeBc=12  #number of basal node
aNodeBc=12 #number of apical nodes
lNodeBc=12

rBc=1

cellCenterPouchX=[domainMinX+cellsGapX+wPouch/2] 
cellCenterPouchY=[domainMinY+hPouchR/2] 
for x in range(nPouchCells-1):
    cellCenterPouchX.append ( cellCenterPouchX[x]+wPouch + cellsGapX)
    cellCenterPouchY.append ( cellCenterPouchY[x]) 


cellCenterPeripX =[domainMinX+cellsGapX+wPerip/2]  
cellCenterPeripY=[domainMinY+hPouchR+lumen+hPerip/2] 
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
    


cellCenterX=[cellCenterPouchX[0]]
cellCenterY=[cellCenterPouchY[0]]

for x in range (1,nPouchCells) : 
    cellCenterX.append (cellCenterPouchX[x] )  
    cellCenterY.append (cellCenterPouchY[x] )  


for x in range (nBcCells) : 
    cellCenterX.append ( cellCenterBcX[x] )  
    cellCenterY.append ( cellCenterBcY[x] )


for x in range (1,nPeripCells+1) :
    cellCenterX.append ( cellCenterPeripX[nPeripCells-x] )  
    cellCenterY.append ( cellCenterPeripY[nPeripCells-x] )  

for x in range (nBcCells,2*nBcCells) : 
    cellCenterX.append ( cellCenterBcX[x] )  
    cellCenterY.append ( cellCenterBcY[x] )


centerPlot(cellCenterX, cellCenterY) 






xC=[ [0 for x in range ( cellTotalNode)] for y in range (len(cellCenterX) ) ]
yC=[ [0 for x in range ( cellTotalNode)] for y in range (len(cellCenterX) ) ] 

for k in range ( nPouchCells ) :
    lastpoint=0 
    if k==0 :
        lNodePouchL=14
        lNodePouchR=42
        hPouchL=2
        hPouchR=6
        dh=2
    elif  k==(nPouchCells-1) :
        lNodePouchL=42
        lNodePouchR=14
        hPouchL=6
        hPouchR=2
        dh=2
    elif k==1 :
        lNodePouchL=42
        lNodePouchR=70
        hPouchL=6
        hPouchR=10
        dh=2
    elif k==(nPouchCells-2) :
        lNodePouchL=70
        lNodePouchR=42
        hPouchL=10
        hPouchR=6
        dh=2
    else:
        lNodePouchL=5*14
        lNodePouchR=5*14
        hPouchL=10
        hPouchR=10
        dh=0

    # from center to top equal to H/2
    for i  in range (int(lNodePouchR/2))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+wPouch/2 
        yC[k][j]=cellCenterY[k]+i/(lNodePouchR/2)*hPouchR/2  

    lastpoint=lastpoint+int (lNodePouchR/2)  
    for i  in range (aNodePouch) :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+wPouch/2- i/aNodePouch*wPouch 
        if k==0 or k==1 :
            yC[k][j]=cellCenterY[k]+hPouchR/2 - dh*(i)/aNodePouch 
        elif k==(nPouchCells-2) or k==(nPouchCells-1) :
            yC[k][j]=cellCenterY[k]+hPouchR/2 + dh*i/aNodePouch  
        else:
            yC[k][j]=cellCenterY[k]+hPouchR/2 

    lastpoint=lastpoint+ aNodePouch  

    for i  in range (lNodePouchL) :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]-wPouch/2 
        yC[k][j]=cellCenterY[k]+hPouchL/2  - i/lNodePouchL*hPouchL 
    
    lastpoint=lastpoint+ lNodePouchL 

    for i  in range (bNodePouch) :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]-wPouch/2+ i/bNodePouch*wPouch

        if k==0 or k==1 :
            yC[k][j]=cellCenterY[k]-hPouchL/2 -dh*i/bNodePouch
        elif k==(nPouchCells-2) or k==(nPouchCells-1) :
            yC[k][j]=cellCenterY[k]-hPouchL/2 +dh*(i)/bNodePouch
        else:
            yC[k][j]=cellCenterY[k]-hPouchL/2 

    lastpoint=lastpoint+ bNodePouch 

    for i  in range (int(lNodePouchR/2))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+wPouch/2 
        yC[k][j]=cellCenterY[k]-hPouchR/2+ i/(lNodePouchR/2)*hPouchR/2  







for k in range ( nPouchCells,nPouchCells+int(nBcCells) ) :
    lastpoint=0
    # from center to top equal to H/2
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/(lNodeBc/2)*np.pi*0.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/(lNodeBc/2)*np.pi*0.25)

    lastpoint=lastpoint+int (lNodeBc/2)  

    for i  in range (int(aNodeBc))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/aNodeBc*np.pi*0.5+np.pi*0.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/aNodeBc*np.pi*0.5+np.pi*0.25)

    lastpoint=lastpoint+int (aNodeBc)  

    for i  in range (int(lNodeBc))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/lNodeBc*np.pi*0.5+np.pi*0.75)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/lNodeBc*np.pi*0.5+np.pi*0.75)

    lastpoint=lastpoint+int (lNodeBc)  

    for i  in range (int(bNodeBc))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/bNodeBc*np.pi*0.5+np.pi*1.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/bNodeBc*np.pi*0.5+np.pi*1.25)

    lastpoint=lastpoint+int (bNodeBc)  

    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/(lNodeBc/2)*np.pi*0.25+np.pi*1.75)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/(lNodeBc/2)*np.pi*0.25+np.pi*1.75)


for k in range ( nPouchCells+int(nBcCells),nPouchCells+int(nBcCells)+nPeripCells ) :

    lastpoint=0
    # from center to top equal to H/2
    for i  in range (int(lNodePerip/2))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+wPerip/2 
        yC[k][j]=cellCenterY[k]+i/(lNodePerip/2)*hPerip/2  

    lastpoint=lastpoint+int (lNodePerip/2)  
    for i  in range (aNodePerip) :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+wPerip/2- i/aNodePerip*wPerip
        yC[k][j]=cellCenterY[k]+hPerip/2  

    lastpoint=lastpoint+ aNodePerip  

    for i  in range (lNodePerip) :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]-wPerip/2 
        yC[k][j]=cellCenterY[k]+hPerip/2  - i/lNodePerip*hPerip 
    
    lastpoint=lastpoint+ lNodePerip 

    for i  in range (bNodePerip) :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]-wPerip/2+ i/bNodePerip*wPerip
        yC[k][j]=cellCenterY[k]-hPerip/2  

    lastpoint=lastpoint+ bNodePerip 

    for i  in range (int(lNodePerip/2))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+wPerip/2 
        yC[k][j]=cellCenterY[k]-hPerip/2+ i/(lNodePerip/2)*hPerip/2  


for k in range ( nPouchCells+nPeripCells+int(nBcCells),nPouchCells+nPeripCells+int(2*nBcCells) ) :
    lastpoint=0
    # from center to top equal to H/2
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/(lNodeBc/2)*np.pi*0.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/(lNodeBc/2)*np.pi*0.25)

    lastpoint=lastpoint+int (lNodeBc/2)  

    for i  in range (int(aNodeBc))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/aNodeBc*np.pi*0.5+np.pi*0.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/aNodeBc*np.pi*0.5+np.pi*0.25)

    lastpoint=lastpoint+int (aNodeBc)  

    for i  in range (int(lNodeBc))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/lNodeBc*np.pi*0.5+np.pi*0.75)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/lNodeBc*np.pi*0.5+np.pi*0.75)

    lastpoint=lastpoint+int (lNodeBc)  

    for i  in range (int(bNodeBc))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/bNodeBc*np.pi*0.5+np.pi*1.25)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/bNodeBc*np.pi*0.5+np.pi*1.25)

    lastpoint=lastpoint+int (bNodeBc)  

    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(i/(lNodeBc/2)*np.pi*0.25+np.pi*1.75)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(i/(lNodeBc/2)*np.pi*0.25+np.pi*1.75)

for k in range ( len(cellCenterX) ) :   
    plt.scatter(xC[k], yC[k])

#plt.xlabel('X (micro meter)')
#plt.ylabel('Y (micro meter)')
#plt.title('locations of cell cetner')
#plt.grid(True)
#plt.savefig("test.png")
plt.show() 

#centerPlot(xC[0], yC[0]) 




