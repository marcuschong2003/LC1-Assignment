# For finding data in gross table
# 1 for Gross Premium
# 2 for Net Premium
# 3 for E(PVFP)
# 4 for E(PVFB)
# 5 for E(PVFE)

import math

def isnumber(a):
    return (isinstance(a,int) or isinstance(a,float))

def Makeham(x,A,B,C):
    return math.exp(-A-B/math.log(C)*math.pow(C,x)*(C-1))

def MakehamMort(A,B,C,Start,End):
    if Start>=End:
        raise Exception("The Starting Age must be strictly smaller than the Ending Age")
    Mortality = []
    for x in range(math.floor(Start),math.ceil(End)+1):
        Mortality.append([x,(A+B*math.pow(C,x)),Makeham(x,A,B,C),1-Makeham(x,A,B,C)])
    return Mortality

def sumproduct(*args):
    product = 0
    length = len(args[0])
    for x in args:
        if len(x)!= length:
            raise Exception("The lenght of all list should be the same")
    for x in args:
        for i in x:
            if not isnumber(i):
                raise Exception("The elements should all be integer or float")
    for i in range(length):
        product += math.prod([args[j][i] for j in range(len(args))])
    return product

class MakehamMortalityData:
    def __init__(self, A,B,C,Start,End):
        self.Start = Start
        self.End = End
        self.Setting = {"A":A,"B":B,"C":C}
        self.data = MakehamMort(A,B,C,Start,End)
        self.span = End-Start
    def find(self,x,index):
        if not (isinstance(x,int) and isinstance(index,int)):
            raise Exception("Age and Index should both be integer")
        if index<1 or index>3:
            raise Exception("Index should be 1 for force of mortality, 2 for survival rate and 3 for death rate")
        if x<self.Start or x>self.End:
            raise Exception("Age not in Range")
        for i in range(len(self.data)):
            if self.data[i][0] == x:
                return self.data[i][index]

def EPVFB(mortality,survivalBenefit,deathBenefit,annualinterest,term,entryage):
    if term>mortality.span:
        raise Exception("Insufficient Mortality data")
    if entryage>(mortality.End - term) or entryage<(mortality.Start):
        raise Exception("Entry Age invalid")
    # if survivalBenefit<deathBenefit:
    #     raise Warning("Survival Benefit should always be larger than death Benefit")
    # return list(range(entryage,entryage+term+1))
    qxs=[mortality.find(x,3) for x in range(entryage+1,entryage+term)]
    vs=[math.pow((1+annualinterest),-x) for x in range(2,term+1)]
    tpxs = [math.prod([mortality.find(x,2) for x in range(entryage,entryage+y)]) for y in range(1,term)]
    nyearterm = sumproduct(qxs,vs,tpxs,deathBenefit[1:])+math.pow((1+annualinterest),-1)*mortality.find(entryage,3)*deathBenefit[0]
    tpxsurvival = math.prod([mortality.find(x,2) for x in range(entryage,entryage+term)])
    APVSurvival = survivalBenefit[-1]*tpxsurvival*math.pow((1+annualinterest),-term)
    epvfb = APVSurvival + nyearterm
    return epvfb

def EPVFP(mortality,survivalBenefit,deathBenefit,annualinterest,term,entryage,frequency=1):
    tpxs = [math.prod([mortality.find(x,2) for x in range(entryage,entryage+y)]) for y in range(1,term)]
    vs=[math.pow((1+annualinterest),-x) for x in range(1,term)]
    d=1-1/(1+annualinterest)
    im=frequency*(math.pow((1+annualinterest),1/frequency)-1)
    dm=im/(1+(im/frequency))
    alpha=annualinterest*d/im/dm
    beta=(annualinterest-im)/(im*dm)
    nyeartermlifeannuity = sumproduct(tpxs,vs)+1
    nEx=math.pow((1+annualinterest),-term)*math.prod([mortality.find(x,2)for x in range(entryage,entryage+term)])
    mthlynyeartermlifeannuity = alpha*nyeartermlifeannuity-beta*(1-nEx)
    epvfb = EPVFB(mortality,survivalBenefit,deathBenefit,annualinterest,term,entryage)
    P = epvfb/nyeartermlifeannuity
    pi = epvfb/mthlynyeartermlifeannuity/frequency
    return {"AnnualPremium":P,"m-thly Premium":pi,"axn":nyeartermlifeannuity,"axnm":mthlynyeartermlifeannuity}

def EPVFE(mortality,survivalBenefit,deathBenefit,annualinterest,term,entryage,frequency=12):
    qxs=[mortality.find(x,3) for x in range(entryage+1,entryage+term)]
    vs=[math.pow((1+annualinterest),-x) for x in range(2,term+1)]
    tpxs = [math.prod([mortality.find(x,2) for x in range(entryage,entryage+y)]) for y in range(1,term)]
    nyearterm = sumproduct(qxs,vs,tpxs)+math.pow((1+annualinterest),-1)*mortality.find(entryage,3)
    epvfb = EPVFB(mortality,survivalBenefit,deathBenefit,annualinterest,term,entryage)
    epvfp = EPVFP(mortality,survivalBenefit,deathBenefit,annualinterest,term,entryage,frequency)
    axn = epvfp.get("axn")
    axnm = epvfp.get("axnm")
    newvs = [math.pow((1+annualinterest),-1)]+vs[:-1]
    expenses = [50]+[4*math.pow(1.02,x)for x in range(term-1)]
    APVperpolicy = sumproduct(newvs,expenses[1:],tpxs)+expenses[0]
    settlementAPV = 100*nyearterm
    GrossPremium = (epvfb + APVperpolicy + settlementAPV)/(axnm-(0.05*axn+(0.15)))
    epvfpamount = GrossPremium * axnm
    epvfeamount = epvfpamount-epvfb
    return {"GrossPremium":GrossPremium,"NetPremium":epvfp.get("AnnualPremium"),"EPVFP":epvfpamount,"EPVFB":epvfb,"EPVFE":epvfeamount}

class Model:
    def __init__(self,A,B,C,StartAge,EndAge,survivalBenefit,deathBenefit,annualinterest,term,frequency):
        self.mortality = MakehamMortalityData(A,B,C,StartAge,EndAge+22)
        self.range = [StartAge,EndAge]
        self.SB = survivalBenefit
        self.DB = deathBenefit
        self.interest = annualinterest
        self.term = term
        self.frequency = frequency
        self.grosstable = [[x,EPVFE(self.mortality,self.SB,self.DB,self.interest,self.term,x,12).get("GrossPremium"),EPVFP(self.mortality,self.SB,self.DB,self.interest,self.term,x,12).get("AnnualPremium"),EPVFE(self.mortality,self.SB,self.DB,self.interest,self.term,x,12).get("EPVFP"),EPVFE(self.mortality,self.SB,self.DB,self.interest,self.term,x,12).get("EPVFB"),EPVFE(self.mortality,self.SB,self.DB,self.interest,self.term,x,12).get("EPVFE")]for x in range(self.range[0],self.range[1]+1)]
    def findage(self,x):
        if not (isinstance(x,int)):
            raise Exception("Age and Index should be number")
        for i in range(len(self.grosstable)):
            if(x==self.grosstable[i][0]):
                return self.grosstable[i]
    def findspecific(self,x,index):
        age = self.findage(x)
        if index<1 or index>5:
            raise Exception("Index should be 1 for Gross Premium, 2 for Net Premium, 3 for E(PVFP), 4 for E(PVFB), 5 for E(PVFE)")
        return age[index]
    def Grosstable(self):
        print(f"{'AGE':<3}|{'Gross Premium':<15}|{'NET Premium':<15}|{'E(PVFP)':<15}|{'E(PVFB)':<15}|{'E(PVFE)':<15}")
        for x in self.grosstable:
            print(f"{x[0]:0>3}|{x[1]:->12}|{x[2]:->12}|{x[3]:->12}|{x[4]:->12}|{x[5]:->12}")
answer = Model(0.03,0.0002,1.06,40,100,[130000]*20,[90000]*20,0.05,20,12)
answer.Grosstable()
