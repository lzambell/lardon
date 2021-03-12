import os,sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir) 
import channelmapper as cm


#example
# IN (CRP, VIEW, CHANNEL)
myinput = [(0, 1, x) for x in range(0, 10)]

print("in (CRP, view, ch) format: ")
print(myinput)

inDAQ   = [cm.CRPToDAQ(crp, v, ch) for crp, v, ch in myinput]
print("\nin DAQ format: ")
print(inDAQ)
