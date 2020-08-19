
from sys import stderr
from contextlib import contextmanager

class LoadingBar:
    
    def __init__(self, label, increments, length = 25, cold_start = False):
        self.increments = increments
        self.length = length
        self.label = label
        self.progress = 0
        self.cold_start = cold_start
        
    def __str__(self):
        if self.cold_start:
            self.cold_start = False
        else:
            self.increment()
        completed_steps = int(self.progress / self.increments * self.length)
        if completed_steps >= self.length:
            return '{}: [{}]'.format(self.label, "="*completed_steps) + '\n' if self.is_finished() else ''
        else:
            return '{}: [{}>{}]'.format(self.label, "="*completed_steps, " "*(self.length - completed_steps - 1))
    
    def increment(self):
        if not self.is_finished():
            self.progress += 1
        
    def is_finished(self):
        return self.progress >= self.increments


class Log:

    def __init__(self, target = stderr):
        self.target = target
        self.indents = 0

    @contextmanager
    def section(self, header):
        try:
            self.start_section(header)
            yield self
        finally:
            self.end_section()

    def start_section(self, section_header):
        self.append(section_header)
        self.indents += 1

    def end_section(self):
        self.indents -= 1 

    def append(self, text, end = '\n', update_line = False):
        linestart = '\r' if update_line else ''
        print(linestart + '\t'*self.indents + str(text), 
            end = '' if update_line else end, 
            file = self.target)