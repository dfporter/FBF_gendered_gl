import re, importlib
import simple.targets

importlib.reload(simple.targets)


class configuration():
    
    def __init__(self, config_text):
        self.settings = {}
        lines = [x.rstrip('\n') for x in config_text.split('\n')]
        in_sublist = False
        for line in lines:
            
            line = line.split('#')[0]
            
            # Do we enter a sublist?
            m_sub = re.search('([^=]+)=\s*\[<(.*)', line)
            if m_sub is not None:
                self.settings[m_sub.group(1)] = {}
                in_sublist = m_sub.group(1)
                continue
            
            # Do we exit a sublist?
            m_end_sub = re.search('(.*)>\]', line)
            if m_end_sub is not None:
                in_sublist = False
                continue
            
            # Do we set a value?
            eq = line.split('=')
            if (len(eq) > 1) and not in_sublist:
                self.settings[re.sub('\s+', '', eq[0])] = '='.join(eq[1:])
            if (len(eq) > 1) and in_sublist:
                self.settings[in_sublist][re.sub('\s+', '', eq[0])] = '='.join(eq[1:])
        
        self.s = self.settings
        
    def targets_from_config(self, name):
        name = str(name)
        t = simple.targets.targets()
        t.read_input(self.s[name])
        return t
