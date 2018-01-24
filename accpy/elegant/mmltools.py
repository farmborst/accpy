from __future__ import print_function,division
import numpy as np
import sys
from scipy import io

class Namemap():
    # A global table with
    # AT index (python convention) - AT ('familiy')name - AO name - PS name - ATtype (defined in AO!)
    def __init__(self,N):
        # NOTE: starts at ZERO (python style), not ONE (matlab style)
        self.ATindex = np.zeros(N, dtype=np.int32)
        # Note, AOindex is not unique but is defined in each familiy (for each ATname)
        self.AOindex = np.zeros(N, dtype=np.int32)
        self.ATtype = np.empty([N], dtype="S40")
        self.ATname = np.empty([N], dtype="S40")
        self.AOname = np.empty([N], dtype="S40")
        self.PSname = np.empty([N], dtype="S40")
        for i in range(N):
            self.ATtype[i] = ''
            self.ATname[i] = ''
            self.AOname[i] = ''
            self.PSname[i] = ''
            
    def getATindicesByAOname(self,name):
        return self.ATindex[self.AOname == name]

    def getATindicesByPSname(self,name):
        return self.ATindex[self.PSname == name]

    def getATindicesByATname(self,name):
        return self.ATindex[self.ATname == name]

    def getPSnames(self,ATtype = 'QUAD'):
        sel = self.ATtype == ATtype
        return np.unique(self.PSname[sel])
    
    def getAOnames(self,ATtype = 'QUAD'):
        sel = self.ATtype == ATtype
        return np.unique(self.AOname[sel])

    def getATnames(self,ATtype = 'QUAD'):
        sel = self.ATtype == ATtype
        return np.unique(self.ATname[sel])

class PrintATRing():
    # can deal with a matlab AT RING objection
    def __init__(self,r, details='default'):
        NATelements = len(r)
        #print(r.shape) #, r[4][0][0][0]
        # (1241,) ([u'BPM'], [[0]], [u'IdentityPass'], [[1700000000]])
        # def getstrings(l):
        #     st = []
        #     for e in r[i][0][0]:
        #         if isinstance(e[0], basestring):
        #             st.append(e[0])
        #     return st
        s = 0
        slast = 0
        for i in range(NATelements):
            # st = getstrings(r[i][0][0])
            # name = st[0]
            # typename = st[1]
            name = r[i][0,0].FamName[0]
            typename = r[i][0,0].PassMethod[0]
            l = r[i][0,0].Length[0,0]
            s = s+l
            #e = r[i][0][0][3][0][0]
            #    print i, name, s, typename,  e,
            # if typename == 'DriftPass':
            #     continue
            # if name in ('VCM','HCM','BPM','HKCM'):
            #     continue
            # if l == 0:
            #     continue
            if i >= 0: # s <= 320 and 
                extra = ''
                if typename in ('StrMPoleSymplectic4Pass','StrMPoleSymplectic4RadPass'):
                    #K = r[i][0][0][2][0][0]
                    try:
                        K = r[i][0,0].K[0,0]
                        extra+=' K = {0:11.8f}'.format(K)
                    except:
                        K = r[i][0,0].PolynomB[0,2]
                        extra+=' PolynomB[0,2] = {0:11.8f}'.format(K)
                if typename in ('ThinMPolePass'):
                    K = r[i][0,0].PolynomB[0,2]
                    extra+=' PolynomB[0,2] = {0:11.8f}'.format(K)
                if typename in ('BndMPoleSymplectic4Pass','BndMPoleSymplectic4RadPass'):
                    # BendingAngle = r[i][0][0][4][0][0]
                    # EntranceAngle = r[i][0][0][5][0][0]
                    # ExitAngle = r[i][0][0][6][0][0]
                    BendingAngle = r[i][0,0].BendingAngle[0,0]
                    EntranceAngle = r[i][0,0].EntranceAngle[0,0]
                    ExitAngle = r[i][0,0].ExitAngle[0,0]
                    extra+=' BendingAngle = {0:11.8f} EntranceAngle = {1:11.8f} ExitAngle = {2:11.8f}'.format(BendingAngle, EntranceAngle, ExitAngle)

                    
                if extra == '' and details == 'default':
                    continue
                print('start: {0:12.6f} - {1:10.6f}   length: {2:10.5f}   '.format(s-l,s,l), '{0:5n} {1:12}  {2:20}'.format(i,name,typename), extra)
                #getstrings(r[i][0][0])
                # r[i][0][0], 
               #, typename 
                #, ' center last end to this start: {0:12.6f}'.format((s-l+slast)*0.5)
                #r[i]
                slast = s

        print('Total length: {0:12.6f}'.format(s))




class ATRingWithAO():
    d = None
    ao = None
    ad = None
    rings = None
    namemap = None
    NATelements = None
    
    def __init__(self, filename):
        # struct_as_record=False preserves nested dictionaries!
        self.d = io.loadmat(filename, struct_as_record=False, squeeze_me=False)

        self.ao = self.d['ao'][0][0]
        self.ad = self.d['ad'][0][0]
        self.rings = self.d['RINGs'][0,:]
        r = self.rings[0].ring[0,:]

        self.NATelements = len(r)
        self.namemap = Namemap(self.NATelements)  

        # loop over AO entries to fill namemap
        for field in self.ao._fieldnames:
            # getattr will return the value for the given key
            #print 'FAM field',field
            if field in ('TUNE','DCCT'):
                continue
            fam = getattr(self.ao, field)[0,0]
            # print fam.Monitor.shape
            # continue
            ATindices = fam.AT[0,0].ATIndex
            ATtype = fam.AT[0,0].ATType[0]
            AOnames = fam.CommonNames
            PSnames = fam.Monitor[0,0].ChannelNames
            Nfirst, Nsecond = ATindices.shape
            #print ATindices.shape, AOnames.shape, PSnames.shape
            # if field != 'BEND':
            #     continue
            # print field, len(ATindices), len(AOnames), len(PSnames)
            # if field == 'BEND':
            #     print ATindices, AOnames, PSnames, ATindices[0] - 1 , ATindices[0] 
            if field == 'RF':
                AOnames = np.repeat(AOnames, len(ATindices))
                PSnames = np.repeat(PSnames, len(ATindices))
                #print ATindices, AOnames, PSnames
            for isecond in range(Nsecond):
                for i in range(Nfirst):
                    j = ATindices[i,isecond] - 1 
                    self.namemap.ATindex[j] = j # note: start at zero (python style)
                    self.namemap.AOindex[j] = i # note: start at zero (python style)
                    self.namemap.ATtype[j] = ATtype
                    self.namemap.AOname[j] = AOnames[i]
                    self.namemap.PSname[j] = PSnames[i].split(':')[0]
                    #print isecond, j,  ATindices[i],  AOnames[i],  PSnames[i]

        # loop over RING to fill additionally with AT ('familiy') name    
        for i in range(self.NATelements):
             self.namemap.ATname[i] = r[i][0,0].FamName[0]

        # some consistency checks
        broken = 0
        for i in range(self.NATelements):
            if self.namemap.ATname[i] == 'BEND':
                if self.namemap.PSname[i] not in ('PB1ID6R','PB2ID6R','PB3ID6R', 'BPR','BPRP'):
                    print(self.namemap.PSname[i])
                    broken = 1
        if broken:
            print('ERROR : BROKEN FILE!!! Probable cause: init and AT file incompatible!!!')
            sys.exit()

    def printNamemap(self, Nmax=None):
        if Nmax == None:
            Nmax = self.NATelements
        nm = self.namemap
        print(' {0:5} {1:5} {2:5}  {3:12}    {4:12}   {5:20}   {6:20}'.format('index','ATindex','AOindex','ATname', 'AOname', 'PSname','ATtype (defined in AO)'))
        for i in range(Nmax):
            print('{0:6n}  {1:6n}  {2:6n}  {3:12}    {4:12}   {5:20}   {6:20}'.format(i,nm.ATindex[i],nm.AOindex[i],nm.ATname[i],nm.AOname[i],nm.PSname[i],nm.ATtype[i]))

    def getRing(self, fitIteration='last'):
        if fitIteration == 'last':
            fitIteration = len(self.rings)-1
        return self.rings[fitIteration].ring[0,:]


    def _printMagnetStrength(self, name, strength, ATtype, outputstyle):
        # print(name)
        # print('Locofile belongs to this machine =', self.ad.Machine)
        Qlength = {'Q1' : 0.25 ,'Q2' : 0.20 ,'Q3' : 0.25 ,'Q4' : 0.50,'Q5' : 0.20}
        Slength = {'S1' : 0.21 ,'S2' : 0.16 ,'S3' : 0.16 ,'S4' : 0.16}
        if self.ad.Machine == 'MLS':
					Qlength = {'Q1' : 0.2 ,'Q2' : 0.2 ,'Q3' : 0.2}
					Slength = {'S1' : 0.1 ,'S2' : 0.1 ,'S3' : 0.1}
        if ATtype == 'SEXT':
            # n = name.replace('PR', 'PR').replace('PD', 'PD').replace('PT', 'PT').replace('PQ', 'PQ').replace('RP', 'RP')
            # n = name.replace('PR', '').replace('PD', 'D').replace('PT', 'T').replace('PQ', 'Q').replace('R', '')   #BESSY
            # n = name.replace('RP', '')			#MLS
            if self.ad.Machine == 'BESSYII':
							n = name.replace('PR', '').replace('PD', 'D').replace('PT', 'T').replace('PQ', 'Q').replace('R', '')
            if self.ad.Machine == 'MLS':
							n = name.replace('RP', '')

            l = Slength[name[:2]]
            strength = strength / l * 2
            if outputstyle=='elegant':
                if n == 'S1':
                     print('! Note: Is S1 split? -> change length to 0.105')
                print('{0:8}: KSEXT, n_kicks="sextkicks", L={1:<8}, K2={2:11.8f}'.format(n,l,strength))
            elif outputstyle=='MADX':
                if n == 'S1':
                     print('// Note: Is S1 split? -> change length to 0.105')
                print('{0:8}: SEXTUPOLE, L={1:<8}, K2={2:11.8f};'.format(n,l,strength))
            elif outputstyle=='visual':
                print('{0:12} K2 = {1:11.8f}'.format(name,strength))
            else:
                print('Output style not implemented.')
        if ATtype == 'QUAD':
            # n = name.replace('PD', 'PD').replace('PT', 'PT').replace('PQ', 'PQ').replace('RP', 'RP')
            # n = name.replace('PD', 'D').replace('PT', 'T').replace('PQ', 'Q').replace('R', '')		#BESSY
            # n = name.replace('RP', '')	#MLS
            if self.ad.Machine == 'BESSYII':
							n = name.replace('PD', 'D').replace('PT', 'T').replace('PQ', 'Q').replace('R', '')
            if self.ad.Machine == 'MLS':
							n = name.replace('RP', '')

            if name[:1] == 'Q':
                l = Qlength[name[:2]]
            if name == 'PQIT6R':
                l = 0.12200 # half length!

            if outputstyle=='elegant':
                print('{0:8}: KQUAD, n_kicks="quadkicks", L={1:<8}, K1={2:11.8f}'.format(n,l,strength))
            elif outputstyle=='MADX':
                print('{0:8}: QUADrupole, L={1:<8}, K1={2:11.8f};'.format(n,l,strength))
            elif outputstyle=='visual':
                print('{0:12} K1 = {1:11.8f}'.format(name,strength))
            else:
                print('Output style not implemented.')
                
    def getMagnetStrength(self, ATtype='QUAD', fitIteration='last', method='byPowerSupply', outputstyle='visual'):
        
        print('Locofile belongs to this machine =', self.ad.Machine)
        if fitIteration == 'last':
            fitIteration = len(self.rings)-1
        r = self.rings[fitIteration].ring[0,:]

        if ATtype == 'QUAD':
            def getS(i):
                return r[i][0,0].K[0,0]
        elif ATtype == 'SEXT':
            def getS(i):
                return r[i][0,0].PolynomB[0,2]
        else:
            print('ATtype not implemented.')
            
        if  method == 'byPowerSupply':
            print('List magnet ({0}) strength by power supply.'.format(ATtype))
            PSn = self.namemap.getPSnames(ATtype = ATtype)
            print('Number of independent parameters:',len(PSn))
            for n in PSn:
                ATi = self.namemap.getATindicesByPSname(n)
                S = getS(ATi[0])
                Savg = 0
                for i in ATi:
                    if S != getS(i): #r[i][0,0].K[0,0]:
                        print('WARNING: Differnt elements of same power supply have differnt values! Probably not fitted according to Power supply! Using average value!')
                    Savg += getS(i)
                Savg = Savg / len(ATi)
                self._printMagnetStrength(n,Savg,ATtype,outputstyle)
        else:
            print('Method not implemented.')


    def getArchiverScalar(self,var,t, index):
        # from AS archiver.py module and modified
        from urllib2 import urlopen,quote
        import datetime

        def filter_camonitor(data):
            while True:
                a,=next(data),
                if not a.startswith('#'):
                    s=a.split()
                    if len(s)==3:
                        s+=['0']
                    yield ' '.join([s[0],s[1]+'+'+s[2],s[3]])
                else:
                    yield a

        def archiver(var,t0,t1,archive='master'):
            #bii="http://archiver.bessy.de/archive/cgi/CGIExport.cgi?INDEX=/opt/Archive/master_index&COMMAND=Raw+Data"
            if archive == 'current_week':
                bii="http://archiver.bessy.de/archive/cgi/CGIExport.cgi?INDEX=/opt/Archive/current_week/index&COMMAND=camonitor"
            else:
                if archive != 'master' :
                    print('WARNING unknown archive',archive)
                bii='http://archiver.bessy.de/archive/cgi/CGIExport.cgi?INDEX=/opt/Archive/' + index + '&COMMAND=camonitor'
        
            if self.ad.Machine == 'MLS':
                bii="http://arc31c.trs.bessy.de/MLS/cgi/CGIExport.cgi?INDEX=%2Fopt%2FArchive%2Fmaster_index&COMMAND=camonitor"
                
            if type(var) is list:
                var='\n'.join(var)

            var="&NAMES=%s"%quote(var)
            spec="&STRSTART=1&STARTSTR=%s&STREND=1&ENDSTR=%s"%(quote(t0),quote(t1))
            res=urlopen(bii+spec+var)
            #print bii+spec+var

            #return np.genfromtxt(filter_Raw(res),dtype=None,skip_footer=1,names=('time','value'),converters={0: lambda d : datetime.datetime.strptime(d[:-3],'%d/%m/%Y+%H:%M:%S.%f')})

            return np.genfromtxt(filter_camonitor(res),dtype=None,names=('name','time','value'),converters={1: lambda d : datetime.datetime.strptime(d[:-3],'%Y-%m-%d+%H:%M:%S.%f')})

        def archiverScalar(var,t):
            return float(archiver(var,t,t)['value'])

        return archiverScalar(var,t)
        
    
    def getMagnetStrengthOnline(self, source='archiver', time='2016-07-28 11:00:00', ATtype='QUAD', method='byPowerSupply', outputstyle='visual', index='master_index'):
        if self.ad.Machine == 'MLS':
            suffix = ':setCur'
            suffixstat1 = ':stPower'
        else:
            suffix = ':set'
            suffixstat1 = ':stat1'
        if source == 'epics':
            from epics import PV
            import datetime
            time = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

        if  method == 'byPowerSupply':
            print('List magnet ({0}) strength by power supply. Source: {1} (time: {2}).'.format(ATtype, source,time))
                
            PSn = self.namemap.getPSnames(ATtype = ATtype)
            print('Number of independent parameters:',len(PSn))
            for n in PSn:
                # take first element
                field = self.namemap.ATname[self.namemap.getATindicesByPSname(n)][0]
                conversionfactors = getattr(self.ao, field)[0,0].Monitor[0,0].HW2PhysicsParams[:,0]
                if self.ad.Machine == 'MLS':
                    # TODO: not understood why needed!?
                    conversionfactors = conversionfactors[0][:,0]   
                # print(conversionfactors)
                cfac = np.mean(conversionfactors)
                if not np.array_equal(conversionfactors, conversionfactors[0]*np.ones_like(conversionfactors)):
                    print('Warning: Different conversion factors for a single power supply given! Taking average.')

                if source=='archiver':
                    Iarch = self.getArchiverScalar(n+suffix, time, index)
                    Savg = cfac * Iarch
                    stat1 = self.getArchiverScalar(n+suffixstat1, time, index)
                    print('cfac*Iarch*stat/L*2 = {}*{}*{}/L*2 = {}/L'.format(cfac, Iarch, stat1, cfac*Iarch*stat1*2))
                    Savg *= stat1
                elif source=='epics':
                    Savg = cfac * PV(n+suffix).get()
                else:
                    print('Source not implemented.')
                    
                self._printMagnetStrength(n,Savg,ATtype,outputstyle)
        else:
            print('Method not implemented.')


