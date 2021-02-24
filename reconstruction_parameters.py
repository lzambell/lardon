import yaml

import config as cf

class Reco:
    def __init__(self, fparam):    
        with open(cf.default_reco) as f: 
            self.param = yaml.load(f, Loader=yaml.FullLoader)

            if(fparam != cf.default_reco):
                with open(fparam) as f: 
                    custom = yaml.load(f, Loader=yaml.FullLoader)
                
                for key, val in custom.items():
                    self.param[key] = val 



