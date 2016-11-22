from ROOT import *
from array import array

import sys


class TreeWriter:

    def __init__(self, varns, N = 0):
        self.varns = varns
        self.N = N
        self.t = None
        self.f = None
        self.vd = [{}, {}]
        self.nvars = None
        self.diffs = None
        
        for vn in self.varns:
            if vn == "event":
                self.vd[0][vn] = array('i', [0])
                self.vd[1][vn] = array('i', [0])
            else:
                self.vd[0][vn] = array('d', [0.0])
                self.vd[1][vn] = array('d', [0.0])
        self.f = TFile.Open("f_"+str(self.N)+".root", "RECREATE")
        self.t = TTree("t", "t")
        for key in self.varns:
            self.t.Branch(key+"_0", self.vd[0][key], key+"_0"+"/"+str(self.vd[0][key].typecode).upper())
            self.t.Branch(key+"_1", self.vd[1][key], key+"_1"+"/"+str(self.vd[1][key].typecode).upper())
        self.nvars = array('i', [len(self.varns)])
        self.t.Branch("nvars", self.nvars, "nvars/I")
        self.diffs = array('i', [0]*len(self.varns))
        self.t.Branch("diffs", self.diffs, "diffs[nvars]/I")

    def fillTree(self, data):
        for ndiff, key in enumerate(self.varns):
            self.vd[0][key][0] = 0
            self.vd[1][key][0] = 0
            self.diffs[ndiff] = 0
        for d in data:
            #print (d[0][2], d[1][2])
            #print self.vd[1]
            #print "Setting the value:", d, self.vd
            self.vd[0][d[0][0]][0] = d[0][2]
            self.vd[1][d[1][0]][0] = d[1][2]
            #print "Getting the value:", d, self.vd
            #print self.vd[1]
            #locals()[d[0][0]+"_0"]
            self.diffs[self.varns.index(d[0][0])] = 1
        if sum(self.diffs) > 0:
        #    print
        #    print self.diffs
        #    print "CIRKOVIC: Filling"
        #    print data
        #    print self.vd[0][key]
        #    print self.vd[1][key]
        #    print
            self.t.Fill()

    def writeFile(self):
        self.t.Write()
        self.f.Close()
'''
def checkIfMerged(vals):
    if True:
        if vals[0].count('.') > 0:
            v = vals[0]
            vals[0] = v[:6]
            vals.insert(1, v[6:])
            #print vals
    #if True:
    #    for j,v in enumerate(vals):
    #        if v.count('.') == 2:
    #            if vals[j][0] != '-':
    #                vals[j] = v[:7]
    #                vals.insert(j+1, v[7:])
    #            else:
    #                vals[j] = v[:8]
    #                vals.insert(j+1, v[8:])
    if True:
        for j,v in enumerate(vals):
            if v.count('.') == 2:
                #print v
                #print v.find('.')
                #print v[:v.find('.')+5]
                #print v[v.find('.')+6:]
                vals[j] = v[:v.find('.')+5]
                vals.insert(j+1, v[v.find('.')+6:])
'''

class PhysicsObject:
    def __init__(self, l, sf="lepton"):
        self.l = ' '.join(l.split())
        #print self.l
        vals = self.l.split(",")
        #checkIfMerged(vals)
        self.props = []
        if sf == "lepton":
            #self.props.append(("evt", "%10d", int(vals[len(self.props)])))
            #self.props.append(("pt", "%10.5f", float(vals[len(self.props)])))
            #self.props.append(("eta", "%10.5f", float(vals[len(self.props)])))
            #self.props.append(("phi", "%10.5f", float(vals[len(self.props)])))
            #self.props.append(("E", "%10.5f", float(vals[len(self.props)])))
            self.props.append(("run", "%6d", int(vals[len(self.props)])))
            self.props.append(("lumi", " %6d", int(vals[len(self.props)])))
            self.props.append(("id", " %10d", int(vals[len(self.props)])))
        #elif sf == "jet":
        #    self.props.append(("evt", "%6d", int(vals[len(self.props)])))
        #    self.props.append(("pt", "%10.5f", float(vals[len(self.props)])))
        #    self.props.append(("eta", "%10.5f", float(vals[len(self.props)])))
        #    self.props.append(("phi", "%10.5f", float(vals[len(self.props)])))
        #    self.props.append(("E", "%10.5f", float(vals[len(self.props)])))
        #elif sf == "tau":
        #    self.props.append(("evt", "%6d", int(vals[len(self.props)])))
        #    self.props.append(("pt", "%10.5f", float(vals[len(self.props)])))
        #    self.props.append(("eta", "%10.5f", float(vals[len(self.props)])))
        #    self.props.append(("phi", "%10.5f", float(vals[len(self.props)])))
        #    self.props.append(("E", "%10.5f", float(vals[len(self.props)])))

    def getPrintString(self):
        '''
        %10d%10.5f%10.5f%10.5f%10.5f%5d%5d%15.5f%15.5f%15.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%21.5f
        '''
        s = ""
        for p in self.props:
            s += (" "+(p[1] % p[2]))
        return s
    
    def getParameter(self, name):
        ret = None
        for p in self.props:
            if p[0] == name:
                ret = p
                break
        return ret

    def getParVal(self, name):
        ret = None
        for p in self.props:
            if p[0] == name:
                ret = p[2]
                break
        return ret

class Lepton(PhysicsObject):
    def __init__(self, l):
        PhysicsObject.__init__(self, l, sf="lepton")
        vals = self.l.split(",")
        #print vals
        #checkIfMerged(vals)
        #self.props.append(("pdgId", "%5d", int(vals[len(self.props)])))
        #self.props.append(("charge", "%5d", int(vals[len(self.props)])))

        #self.props.append(("miniRelIso", "%15.5f", float(vals[len(self.props)])))
        ##self.props.append(("miniRelIsoCharged", "%15.5f", float(vals[len(self.props)])))
        #self.props.append(("miniIsoCharged", "%15.5f", float(vals[len(self.props)])))
        ##self.props.append(("miniRelIsoNeutral", "%15.5f", float(vals[len(self.props)])))
        #self.props.append(("miniIsoNeutral", "%15.5f", float(vals[len(self.props)])))
        #self.props.append(("jetPtRel", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("jetCSV", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("jetPtRatio", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("sip3d", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("dxy", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("dz", "%10.5f", float(vals[len(self.props)])))

        #self.props.append(("lepId", "  %+2d", int(vals[len(self.props)])))
        #self.props.append(("lepPt", "  %6.2f", float(vals[len(self.props)])))
        #self.props.append(("lepEta", " %+4.2f", float(vals[len(self.props)])))
        #self.props.append(("lepPhi", " %+4.2f", float(vals[len(self.props)])))
        #self.props.append(("metpt", "    %6.1f", float(vals[len(self.props)])))
        #self.props.append(("metphi", "  %+4.2f", float(vals[len(self.props)])))
        #self.props.append(("njets", "    %d", int(vals[len(self.props)])))
        #self.props.append(("nbjets", " %d", int(vals[len(self.props)])))
        '''
        self.props.append(("superClusterEta", " %6.3f", float(vals[len(self.props)])))
        self.props.append(("deltaEtaIn", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("deltaPhiIn", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("sigmaIEtaIEta_full5x5", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("hadronicOverEm", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("ioEmIoP", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("pfElectronIso", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("missingHits", " %6d", int(vals[len(self.props)])))
        '''

class Jet(PhysicsObject):
    def __init__(self, l):
        PhysicsObject.__init__(self, l, sf="jet")
        vals = self.l.split(" ")
        #print vals
        #checkIfMerged(vals)
        #self.props.append(("CSVv2", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("MET pt", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("MET phi", "%10.5f", float(vals[len(self.props)])))

class Tau(PhysicsObject):
    def __init__(self, l):
        PhysicsObject.__init__(self, l, sf="tau")
        vals = self.l.split(" ")
        #print vals
        #checkIfMerged(vals)
        #self.props.append(("dxy", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("dz", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("decayModeFinding", "%10.5f", float(vals[len(self.props)])))
        #self.props.append(("byLooseCombinedIsolationDeltaBetaCorr3Hits", "%10.5f", float(vals[len(self.props)])))

class Electron(Lepton):
    def __init__(self, l):
        Lepton.__init__(self, l)
        vals = self.l.split(",")
        #print vals
        #checkIfMerged(vals)
        #self.props.append(("eleMVA", "%21.5f", float(vals[len(self.props)])))
        self.props.append(("superClusterEta", " %6.3f", float(vals[len(self.props)])))
        self.props.append(("deltaEtaIn", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("deltaPhiIn", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("sigmaIEtaIEta_full5x5", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("hadronicOverEm", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("ioEmIoP", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("pfElectronIso", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("missingHits", " %6d", int(vals[len(self.props)])))

        self.props.append(("sumChargedHadronPt", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("sumNeutralHadronEt", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("sumPhotonEt", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("rho", " %6.5f", float(vals[len(self.props)])))
        self.props.append(("Aeff", " %6.5f", float(vals[len(self.props)])))

class Muon(Lepton):
    def __init__(self, l):
        L3epton.__init__(self, l)
        vals = self.l.split(",")
        print vals
        #checkIfMerged(vals)
        #self.props.append(("segmentCompatibility", "%21.5f", float(vals[len(self.props)])))

def main(argv):
    print argv
    roundingErr = 0.01
    with open(argv[1]) as f1:
        content1 = f1.readlines()
    with open(argv[2]) as f2:
        content2 = f2.readlines()
    objs1 = []
    objs2 = []
    if argv[0] == "electrons":
        for l in content1[1:]:
            objs1.append(Electron(l))
            #objs1.append(Lepton(l))
        for l in content2[1:]:
            objs2.append(Electron(l))
            #objs2.append(Lepton(l))
    elif argv[0] == "muons":
        for l in content1[1:]:
            #objs1.append(Muon(l))
            objs1.append(Lepton(l))
        for l in content2[1:]:
            #objs2.append(Muon(l))
            objs2.append(Lepton(l))
    #elif argv[0] == "jets":
    #    for l in content1[1:]:
    #        objs1.append(Jet(l))
    #    for l in content2[1:]:
    #        objs2.append(Jet(l))
    #elif argv[0] == "taus":
    #    for l in content1[1:]:
    #        objs1.append(Tau(l))
    #    for l in content2[1:]:
    #        objs2.append(Tau(l))
    #evts1 = set([o.getParVal("evt") for o in objs1])
    #evts2 = set([o.getParVal("evt") for o in objs2])
    evts1 = set([o.getParVal("id") for o in objs1])
    evts2 = set([o.getParVal("id") for o in objs2])
    ex1 = sorted(list(evts1-evts2))
    ex2 = sorted(list(evts2-evts1))
    ies = sorted(list(evts1.intersection(evts2)))
    print
    #print "<: ", ex1
    #print ">: ", ex2
    #print "=: ", ies

    if argv[0] == "muons":
        print "MUONS"
#        print "absolute dxy, dz diffs skipped"
#        print "jetPtRel rounding diffs skipped"
#        print "jetPtRatio rounding diffs skipped"
        print "all rounding diffs skipped"
#        print "jetCSV == -10.0 vs. jetCSV == 0.0 skipped"
        #print "skip jetCSV comparing"
        #print "E rounding diffs skipped"
        #print "skip jetPtRel comparing"
#        print "skipp all miniIsoNeutral(0.0, < 0.0) diffs"
#        print "miniIsoNeutral rounding diffs skipped"
        print "<: ", len(ex1), ex1
        print ">: ", len(ex2), ex2
        print "=: ", len(ies), ies
        #sys.exit()

        #tw = TreeWriter(["event", "pT", "Eta", "Phi", "EpdgIDcharge", "miniRelIso", "miniIsoCharged", "miniIsoNeutral", "jetPtRel", "jetCSV", "jetPtRatio", "sip3d", "dxy", "dz", "segmentCompatibility"], "muons")
        tw = TreeWriter(["run", "lumi", "id", "lepId", "lepPt", "lepEta", "lepPhi", "metpt", "metphi", "njets", "nbjets"], "muons")

        iesl = list(ies[:int(argv[3])]) if (int(argv[3]) != -1) else ies
        for i in iesl:
            ps1 = [o for o in objs1 if o.getParVal("id") == i][0]
            ps2 = [o for o in objs2 if o.getParVal("id") == i][0]
    #        if len(ps1) > 1 or len(ps2) > 1:
    #            sys.exit()

            data = []

            dprsn = str(i)+" ["
            for p in xrange(0, len(ps1.props)):
                if ps1.props[p] != ps2.props[p]:
                    if False:
                        continue
#                    elif (ps1.props[p][0] in ["dxy", "dz"]) and abs(ps1.props[p][2]) == abs(ps2.props[p][2]): # absolute dxy, dz diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # jetPtRel rounding diffs skipped
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2]<0.02 if ps1.props[p][2] != 0.0 else False): # jetPtRel rounding diffs skipped
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # jetPtRel rounding diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetPtRatio"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # jetPtRatio rounding diffs skipped
#                        continue
                    elif (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # all rounding diffs skipped
                        continue
#                    elif (ps1.props[p][0] in ["jetCSV"]) and (ps1.props[p][2] == -10.0) and (ps2.props[p][2] == 0.0): # jetCSV == -10.0 vs. jetCSV == 0.0 skipped
#                        continue
                    #elif (ps1.props[p][0] in ["jetCSV"]) and (ps1.props[p][2] == -99.0) and (ps2.props[p][2] == 0.0): # jetCSV == -99.0 vs. jetCSV == 0.0 skipped
                    #    continue
                    #elif (ps1.props[p][0] in ["jetCSV"]): # skip jetCSV comparing for muons
                    #    continue
                    #elif (ps1.props[p][0] in ["jetPtRel"]): # skip jetPtRel comparing
                    #    continue
                    #elif (ps1.props[p][0] in ["E"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # E rounding diffs skipped
                    #    continue
                    #elif (ps1.props[p][0] in ["eleMVA"]):
                    #    print ps1.props[p], ps2.props[p]
                    #    sys.exit()
                    #    continue
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and (ps1.props[p][2] == 0.0 and ps2.props[p][2] < 0): # skipp all miniIsoNeutral(0.0, < 0.0) diffs
#                        continue
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # miniIsoNeutral rounding diffs skipped
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # miniIsoNeutral rounding diffs skipped
#                        continue
                    else:
                        #dprsn += " "+ps1.props[p][0]
                        dprsn += " "+ps1.props[p][0]+"("+((ps1.props[p][1]) % (ps1.props[p][2])).replace(" ", "")+", "+((ps2.props[p][1]) % (ps2.props[p][2])).replace(" ", "")+")"

                        data.append((ps1.props[p], ps2.props[p]))

            dprsn += " ]"
            #print ps1[0].getPrintString()
            #print ps2[0].getPrintString()
            #print
            print dprsn

            tw.fillTree(data)

        tw.writeFile()

    elif argv[0] == "electrons":
        print "ELECTRON VARS"
#        print "absolute dxy, dz diffs skipped"
#        print "jetPtRel rounding diffs skipped"
#        print "jetPtRatio rounding diffs skipped"
#        print "all rounding diffs skipped"
#        print "skip jetCSV comparing"
        #print "E rounding diffs skipped"
#        print "skipp all miniIsoNeutral(0.0, < 0.0) diffs"
#        print "miniIsoNeutral rounding diffs skipped"
        print "<: ", len(ex1), ex1
        print ">: ", len(ex2), ex2
        print "=: ", len(ies), ies
        #sys.exit()
        iesl = list(ies[:int(argv[3])]) if (int(argv[3]) != -1) else ies

        #tw = TreeWriter(["event", "pT", "Eta", "Phi", "EpdgIDcharge", "miniRelIso", "miniIsoCharged", "miniIsoNeutral", "jetPtRel", "jetCSV", "jetPtRatio", "sip3d", "dxy", "dz", "eleMVA"], "electrons")
        #tw = TreeWriter(["run", "lumi", "id", "lepId", "lepPt", "lepEta", "lepPhi", "metpt", "metphi", "njets", "nbjets"], "electrons")
        tw = TreeWriter(["run", "lumi", "id", "superClusterEta", "deltaEtaIn", "deltaPhiIn", "sigmaIEtaIEta_full5x5", "hadronicOverEm", "ioEmIoP", "pfElectronIso", "missingHits", "sumChargedHadronPt", "sumNeutralHadronEt", "sumPhotonEt", "rho", "Aeff"], "electrons")

        for i in iesl:
            ps1 = [o for o in objs1 if o.getParVal("id") == i][0]
            ps2 = [o for o in objs2 if o.getParVal("id") == i][0]
    #        if len(ps1) > 1 or len(ps2) > 1:
    #            sys.exit()

            data = []

            dprsn = str(i)+" ["
            for p in xrange(0, len(ps1.props)):
                if ps1.props[p] != ps2.props[p]:
                    if False:
                        continue
#                    elif (ps1.props[p][0] in ["dxy", "dz"]) and abs(ps1.props[p][2]) == abs(ps2.props[p][2]): # absolute dxy, dz diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # jetPtRel rounding diffs skipped
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2]<0.02 if ps1.props[p][2] != 0.0 else False): # jetPtRel rounding diffs skipped
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # jetPtRel rounding diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetPtRatio"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # jetPtRatio rounding diffs skipped
#                        continue
#                    elif (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # all rounding diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetCSV"]): # skip jetCSV comparing
#                        continue
                    #elif (ps1.props[p][0] in ["eleMVA"]):
                    #    print ps1.props[p], ps2.props[p]
                    #    sys.exit()
                    #    continue
                    #elif (ps1.props[p][0] in ["E"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # E rounding diffs skipped
                    #    continue
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and (ps1.props[p][2] == 0.0 and ps2.props[p][2] < 0): # skipp all miniIsoNeutral(0.0, < 0.0) diffs
#                        continue
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # miniIsoNeutral rounding diffs skipped
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # miniIsoNeutral rounding diffs skipped
#                        continue
                    else:
                        #dprsn += " "+ps1.props[p][0]
                        dprsn += " "+ps1.props[p][0]+"("+((ps1.props[p][1]) % (ps1.props[p][2])).replace(" ", "")+", "+((ps2.props[p][1]) % (ps2.props[p][2])).replace(" ", "")+")"

                        data.append((ps1.props[p], ps2.props[p]))

            dprsn += " ]"
            #print ps1[0].getPrintString()
            #print ps2[0].getPrintString()
            #print
            print dprsn

            tw.fillTree(data)

        tw.writeFile()

    elif argv[0] == "jets":
        print "JETS"
        print "all rounding diffs skipped"
        print "<: ", ex1
        print ">: ", ex2
        print "=: ", ies
        iesl = list(ies[:int(argv[3])]) if (int(argv[3]) != -1) else ies

        tw = TreeWriter(["event", "pT", "Eta", "Phi", "E", "CSVv2", "MET pt", "MET phi"], "jets")

        for i in iesl:
            ps1 = [o for o in objs1 if o.getParVal("evt") == i][0]
            ps2 = [o for o in objs2 if o.getParVal("evt") == i][0]
    #        if len(ps1) > 1 or len(ps2) > 1:
    #            sys.exit()

            data = []

            dprsn = str(i)+" ["
            for p in xrange(0, len(ps1.props)):
                if ps1.props[p] != ps2.props[p]:
                    if False:
                        continue
                    elif (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # all rounding diffs skipped
                        continue
#                    elif (ps1.props[p][0] in ["dxy", "dz"]) and abs(ps1.props[p][2]) == abs(ps2.props[p][2]): # absolute dxy, dz diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # jetPtRel rounding diffs skipped
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2]<0.02 if ps1.props[p][2] != 0.0 else False): # jetPtRel rounding diffs skipped
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # jetPtRel rounding diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetPtRatio"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # jetPtRatio rounding diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetCSV"]): # skip jetCSV comparing
#                        continue
                    #elif (ps1.props[p][0] in ["eleMVA"]):
                    #    print ps1.props[p], ps2.props[p]
                    #    sys.exit()
                    #    continue
                    #elif (ps1.props[p][0] in ["E"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # E rounding diffs skipped
                    #    continue
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and (ps1.props[p][2] == 0.0 and ps2.props[p][2] < 0): # skipp all miniIsoNeutral(0.0, < 0.0) diffs
#                        continue
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # miniIsoNeutral rounding diffs skipped
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # miniIsoNeutral rounding diffs skipped
#                        continue
                    else:
                        #dprsn += " "+ps1.props[p][0]
                        dprsn += " "+ps1.props[p][0]+"("+((ps1.props[p][1]) % (ps1.props[p][2])).replace(" ", "")+", "+((ps2.props[p][1]) % (ps2.props[p][2])).replace(" ", "")+")"

                        data.append((ps1.props[p], ps2.props[p]))

            dprsn += " ]"
            #print ps1[0].getPrintString()
            #print ps2[0].getPrintString()
            #print
            print dprsn

            tw.fillTree(data)

        tw.writeFile()

    elif argv[0] == "taus":
        print "TAUS"
        print "all rounding diffs skipped"
        print "<: ", ex1
        print ">: ", ex2
        print "=: ", ies
        iesl = list(ies[:int(argv[3])]) if (int(argv[3]) != -1) else ies

        tw = TreeWriter(["event", "pT", "Eta", "Phi", "E", "dx", "dz", "decayModeFinding", "byLooseCombinedIsolationDeltaBetaCorr3Hits"], "taus")

        for i in iesl:
            ps1 = [o for o in objs1 if o.getParVal("evt") == i][0]
            ps2 = [o for o in objs2 if o.getParVal("evt") == i][0]
    #        if len(ps1) > 1 or len(ps2) > 1:
    #            sys.exit()

            data = []

            dprsn = str(i)+" ["
            for p in xrange(0, len(ps1.props)):
                if ps1.props[p] != ps2.props[p]:
                    if False:
                        continue
                    elif (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # all rounding diffs skipped
                        continue
#                    elif (ps1.props[p][0] in ["dxy", "dz"]) and abs(ps1.props[p][2]) == abs(ps2.props[p][2]): # absolute dxy, dz diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # jetPtRel rounding diffs skipped
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2]<0.02 if ps1.props[p][2] != 0.0 else False): # jetPtRel rounding diffs skipped
#                    elif (ps1.props[p][0] in ["jetPtRel"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # jetPtRel rounding diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetPtRatio"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # jetPtRatio rounding diffs skipped
#                        continue
#                    elif (ps1.props[p][0] in ["jetCSV"]): # skip jetCSV comparing
#                        continue
                    #elif (ps1.props[p][0] in ["eleMVA"]):
                    #    print ps1.props[p], ps2.props[p]
                    #    sys.exit()
                    #    continue
                    #elif (ps1.props[p][0] in ["E"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # E rounding diffs skipped
                    #    continue
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and (ps1.props[p][2] == 0.0 and ps2.props[p][2] < 0): # skipp all miniIsoNeutral(0.0, < 0.0) diffs
#                        continue
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and round(ps1.props[p][2], 2) == round(ps2.props[p][2], 2): # miniIsoNeutral rounding diffs skipped
#                    elif (ps1.props[p][0] in ["miniIsoNeutral"]) and (abs(ps1.props[p][2]-ps2.props[p][2])/ps1.props[p][2] < roundingErr if ps1.props[p][2] != 0.0 else False): # miniIsoNeutral rounding diffs skipped
#                        continue
                    else:
                        #dprsn += " "+ps1.props[p][0]
                        dprsn += " "+ps1.props[p][0]+"("+((ps1.props[p][1]) % (ps1.props[p][2])).replace(" ", "")+", "+((ps2.props[p][1]) % (ps2.props[p][2])).replace(" ", "")+")"

                        data.append((ps1.props[p], ps2.props[p]))

            dprsn += " ]"
            #print ps1[0].getPrintString()
            #print ps2[0].getPrintString()
            #print
            print dprsn

            tw.fillTree(data)

        tw.writeFile()

    #print list(set(evts1).intersection(evts2))

if __name__ == "__main__":
    main(sys.argv[1:])


