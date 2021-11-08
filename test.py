from matplotlib.pyplot import *
import numpy as np 
from netCDF4 import Dataset
import pandas as pd
#from wrf import getvar, interplevel

def findArgmin(data, value):
    IdxMinflat = np.argmin(abs(data - value))
    idxMin = np.unravel_index(IdxMinflat, data.shape)
    return idxMin

class WRFData:
    def __init__(self, fileDir):
        self.filedir = fileDir
        self.wrfdata = Dataset(fileDir)
        self.xlon = self.getVar('XLONG')#self.wrfdata["XLONG"]
        self.xlat = self.getVar('XLAT')#self.wrfdata["XLAT"]
        self.pdTimes = self.constructLocalTime()

    def delSpinUpTime(self, var):
        """remove the spin up (first 24 hours)"""
        var = var[13:]
        return var

    def getVar(self, varName):
        var = np.array(self.wrfdata[varName])
        var = self.delSpinUpTime(var)
        return var

    def UTCtoLocal(self, UTCTimes):
        LocalTimes = UTCTimes + pd.Timedelta(8, 'hour')
        return LocalTimes


    def constructLocalTime(self):
        oriTimes = self.getVar('Times')
        newTimes = []
        for i in range(len(oriTimes)):
            newTimes.append(pd.to_datetime("".join(char.decode("utf-8") for char in oriTimes[i]), format="%Y-%m-%d_%H:%M:%S"))
        newTimes = pd.DatetimeIndex(newTimes)

        LocalTimes = self.UTCtoLocal(newTimes)
        return LocalTimes

class XitouData:
    def __init__(self, fileDir):
        self.wrfdata = WRFData(fileDir)
        self.times = self.wrfdata.pdTimes
        self.xitouLon = 120.7838246
        self.xitouLat = 23.6759616
        self.idxLon, self.idxLat = self.getXitouIdx()
        self.T2 = self.getVarValue("T2") # K
        self.Q2 = self.getVarValue("Q2") # kg / kg
        self.Psrf = self.getVarValue("PSFC") # hPa
        self.ev = self.Qv2Ev() # hPa
        self.evs = self.getEvs() # hPa
        self.RH = self.getRH() 
    def getXitouIdx(self):
        idxGridLon = findArgmin(self.wrfdata.xlon, self.xitouLon)[2]
        idxGridLat = np.argmin(abs(self.wrfdata.xlat[:, :, idxGridLon] - self.xitouLat))
        return idxGridLon, idxGridLat

    def getVarValue(self, varName):
        VarField = self.wrfdata.getVar(varName)
        VarValue = np.array(VarField)[:, self.idxLat, self.idxLon]
        return VarValue

    def Qv2Ev(self):
        ev = self.Psrf/100 * self.Q2 / (self.Q2 + 0.622) # hPa
        return ev

    def getEvs(self):
        """
        use Goff-Gratch, 1946. 
        input: T (K)
        output: es (hPa)
        """
        T = self.T2
        Tst = 373.15 # boiling point (K)
        ln_es = -7.90298 * (Tst / T - 1) + \
             5.02808 * np.log10(Tst / T) - \
             1.3816e-7 * (10 ** (11.344 * (1 - T / Tst)) -1) + \
             8.1328e-3 * (10 ** (-3.49149 * (Tst / T - 1)) - 1) + \
             np.log10(1013.25) # saturated vapor pressure (hPa)
        es = 10 ** (ln_es)
        return es 

    def getRH(self):
        RH = self.ev / self.evs
        return RH

class ModeCollector:
    def __init__(self):
        self.times = []
        self.T2s = []
        self.Q2s = []
        self.RHs = []

    def collectTimes(self, value):
        self.times.append(value)
    def collectT2(self, value):
        self.T2s.append(value)
    def collectQ2(self, value):
        self.Q2s.append(value)
    def collectRH(self, value):
        self.RHs.append(value)

    def collectData(self, PlaceData):
        self.collectTimes(PlaceData.times)
        self.collectT2(PlaceData.T2)
        self.collectQ2(PlaceData.Q2)
        self.collectRH(PlaceData.RH)

    def squeezeList(self):
        self.times = pd.DatetimeIndex(np.hstack(self.times))
        self.T2s = np.hstack(self.T2s)
        self.Q2s = np.hstack(self.Q2s)
        self.RHs = np.hstack(self.RHs)


class DrawSys:
    def __init__(self, ModeName, Mode):
        self.ModeName = ModeName
        self.Mode = Mode

    def drawT2(self):
        figure(figsize=(20, 8))
        grid(True)
        for i, mode in enumerate(self.ModeName):
            plot(self.Mode[i].times, self.Mode[i].T2s-273, label=self.ModeName[i])
        xticks(self.Mode[0].times[::12], ["{MONTH}-{DAY}-{HOUR}".format(MONTH=x.month, DAY=x.day, HOUR=x.hour) for x in self.Mode[0].times[::12]])
        legend()
        ylim(10, 25)
        ylabel(r"[$\degree C$]")
        title("T2 from {T1} to {T2}".format(T1=self.Mode[0].times[0], T2=self.Mode[0].times[-1]))
        savefig("{MODE}_T2.jpg".format(MODE=self.ModeName[0][0]), dpi=300)

    def drawQ2(self):
        figure(figsize=(20, 8))
        grid(True)
        for i, mode in enumerate(self.ModeName):
            plot(self.Mode[i].times, self.Mode[i].Q2s, label=self.ModeName[i])
        xticks(self.Mode[0].times[::12], ["{MONTH}-{DAY}-{HOUR}".format(MONTH=x.month, DAY=x.day, HOUR=x.hour) for x in self.Mode[0].times[::12]])
        legend()
        ylim(0.005, 0.020)
        ylabel(r"[$kg / kg$]")
        title("Q2 from {T1} to {T2}".format(T1=self.Mode[0].times[0], T2=self.Mode[0].times[-1]))
        savefig("{MODE}_Q2.jpg".format(MODE=self.ModeName[0][0]), dpi=300)

    def drawRH2(self):
        figure(figsize=(20, 8))
        grid(True)
        for i, mode in enumerate(self.ModeName):
            plot(self.Mode[i].times, self.Mode[i].RHs*100, label=self.ModeName[i])
        xticks(self.Mode[0].times[::12], ["{MONTH}-{DAY}-{HOUR}".format(MONTH=x.month, DAY=x.day, HOUR=x.hour) for x in self.Mode[0].times[::12]])
        legend()
        ylim(50, 110)
        ylabel("%")
        title("T2 from {T1} to {T2}".format(T1=self.Mode[0].times[0], T2=self.Mode[0].times[-1]))
        savefig("{MODE}_RH.jpg".format(MODE=self.ModeName[0][0]), dpi=300)

if __name__ == "__main__":

    #dateList = [x for x in range(15, 25)]
    dateList = pd.date_range("20210415T12", periods=9, freq="D")
    NC = ModeCollector()
    NM = ModeCollector()
    NU = ModeCollector()
    WC = ModeCollector()
    WM = ModeCollector()
    WU = ModeCollector()
    NModeName = ['NC', 'NM', 'NU']
    WModeName = ['WC', 'WM', 'WU']
    NMode = [NC, NM, NU]
    WMode = [WC, WM, WU]
    for i, mode in enumerate(NModeName):
        print(mode)
        for j, date in enumerate(dateList):
            wrf_dir = "/home/twsand/fskao/wrfOUT43v1/{MODE}202104{DATE}/wrfout_d04_2021-04-{DATE}_12:00:00".format(MODE=mode, DATE=date.day)
            NMode[i].collectData(XitouData(wrf_dir))
        NMode[i].squeezeList()
    Ndraw = DrawSys(NModeName, NMode)
    Ndraw.drawT2()
    Ndraw.drawQ2()
    Ndraw.drawRH2()

    for i, mode in enumerate(WModeName):
        print(mode)
        for j, date in enumerate(dateList):
            wrf_dir = "/home/twsand/fskao/wrfOUT43v1/{MODE}202104{DATE}/wrfout_d04_2021-04-{DATE}_12:00:00".format(MODE=mode, DATE=date.day)
            WMode[i].collectData(XitouData(wrf_dir))
        WMode[i].squeezeList()
    Wdraw = DrawSys(WModeName, WMode)
    Wdraw.drawT2()
    Wdraw.drawQ2()
    Wdraw.drawRH2()



