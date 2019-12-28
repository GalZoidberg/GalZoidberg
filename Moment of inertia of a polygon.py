def J(x,y,M,premik):
    import matplotlib.pyplot as plt
    import math


    def odstrani_nti_element(n,seznam):
        if n==-1:
            seznam=seznam[:len(seznam)-1]
        else:    
            seznam=seznam[:n]+seznam[n+1:]
        return[seznam]

    def ali_se_sekata(x1,x2,x3,x4,y1,y2,y3,y4):
            x=(x2*x4*(y1-y3)+x1*x4*(-y2+y3)+x1*x3*(y2-y4)+x2*x3*(-y1+y4))/(x4*(y1-y2)+x3*(-y1+y2)+(x1-x2)*(y3-y4))
            y=(x4*(y1-y2)*y3+x1*y2*y3-x3*y1*y4-x1*y2*y4+x3*y2*y4+x2*y1*(-y3+y4))/(x4*(y1-y2)+x3*(-y1+y2)+(x1-x2)*(y3-y4))
            if (x1<x<x2 or x2<x<x1) and (x3<x<x4 or x4<x<x3) and (y1<y<y2 or y2<y<y1) and (y3<y<y4 or y4<y<y3):
                return(1)
            else:
                return(0)        #1,2: od ene daljice; 3,4 pa od druge

    def zamakni(množica,premik):        #zamakne točke na navpičnih, vodoravnih premicah
        stop=0
        while stop==0:
            stop=1
            for n in range(len(množica)):
                if množica[n] in množica[:n]:
                    množica[n]+=premik
                    stop*=0
        return(množica)

    def Jtikotnika(a,b,c,m):
        stranice=[a,b,c]
        a2=min(stranice)
        c2=max(stranice)
        b2=stranice[0]+stranice[1]+stranice[2]-a2-c2
        k=math.sqrt(2*b2**2*c2**2  -2*a2**2*b2**2  -2*a2**2*c2**2     +a2**4  +b2**4  +c2**4)
        return(m*k*(2*b2**2+3*c2**2-2*k)/36/c2**2)

    l0=len(x)
    

    dalx2=[x[0],x[-1]]
    daly2=[y[0],y[-1]]
    x=zamakni(x,premik)
    y=zamakni(y,premik)


    trikotniki=[]
    dalx=x
    daly=y
    n=0   ###vse skup preverjamo za n-1 to oglišče
    while len(x)!=3:
        if (x[n-1]-x[n-2])*(y[n]-y[n-1])>(x[n]-x[n-1])*(y[n-1]-y[n-2]):
            sekata=0
            for m in range(len(x)):
                if ali_se_sekata(x[n-2],x[n],x[m-1],x[m],y[n-2],y[n],y[m-1],y[m])==1:
                    sekata=1
            if sekata==0:
                trikotniki+=[[x[n-2],x[n-1],x[n],y[n-2],y[n-1],y[n]]]
                dalx=dalx+[x[n-2],x[n]]
                daly=daly+[y[n-2],y[n]]
                x=odstrani_nti_element(n-1,x)[0]#
                y=odstrani_nti_element(n-1,y)[0]#
        n=(n+1)%(len(x))

    trikotniki+=[x+y]

    težiščax=[]
    težiščay=[]
    S=[]
    for n in trikotniki:
        težiščax+=[(n[0]+n[1]+n[2])/3]
        težiščay+=[(n[3]+n[4]+n[5])/3]
        S+=[abs(n[0]*(n[4]-n[5])+n[1]*(n[5]-n[3])+n[2]*(n[3]-n[4]))/2]

    vsotaS=0
    for n in S:
        vsotaS+=n
    števcx=števcy=0
    for n in range(len(S)):
        števcx+=težiščax[n]*S[n]
        števcy+=težiščay[n]*S[n]
    skupno_težx=števcx/vsotaS
    skupno_težy=števcy/vsotaS

    J=0
    for n in range(len(trikotniki)):
        a=math.sqrt((trikotniki[n][0]-trikotniki[n][1])**2   +   (trikotniki[n][3]-trikotniki[n][4])**2)
        b=math.sqrt((trikotniki[n][0]-trikotniki[n][2])**2   +   (trikotniki[n][3]-trikotniki[n][5])**2)
        c=math.sqrt((trikotniki[n][1]-trikotniki[n][2])**2   +   (trikotniki[n][4]-trikotniki[n][5])**2)
        J+=M*S[n]/vsotaS*(Jtikotnika(a,b,c,1)+(skupno_težx-težiščax[n])**2+(skupno_težy-težiščay[n])**2)


    print('T(%s,%s)'%(str(skupno_težx),str(skupno_težy)))
    print('J=%s'%(str(J)))
    narisanox=narisanoy=[]
    for n in range(0, l0):
        plt.plot(dalx[n:n+2], daly[n:n+2], 'ro-', color='black')
        narisanox+=[[dalx[n],dalx[n+1]]]
        narisanoy+=[[daly[n],daly[n+1]]]
    plt.plot(dalx2, daly2, 'ro-', color='black')
    for n in range(l0, len(dalx)):
        if dalx[n:n+2] not in narisanox or daly[n:n+2] not in narisanoy: 
            plt.plot(dalx[n:n+2], daly[n:n+2], 'ro-')
    plt.plot(težiščax, težiščay, 'ro',color='green')
    plt.plot(skupno_težx, skupno_težy, 'ro',color='blue')
    plt.show()



#PODATKI
#####################################################################################################
x=[-2,0,-1.5,1.5,5,11,4,-3,0]  #x koordniate oglišč naštetih v pozitivni orientaciji                           
y=[-17,-18,-15,-4,-10,9,12,21,0] #y koordniate oglišč naštetih v pozitivni orientaciji
M=1237            #masa prizme
premik=0.0000001  #mejhn premik da ni deljenja z 0, če sta točki na isti navpični premici
######################################################################################################
J(x,y,M,premik)
