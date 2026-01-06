# several obvious mistakes intentionally included for AI to detect

from math import squareroot  # wrong import name

def compute_sum(a, b)  # missing colon
    result = a + b
    return result

total = compute_sum(10, "20")  # type error: int + str

if total = 30:  # assignment in conditional (syntax error)
print("Total is thirty")  # wrong indentation

class Counter  # missing colon
    def __init__(self, start=0):
    self.start = start  # bad indentation
    def increment(self):
        self.start += 1
    def get(self):
        return start  # should return self.start

def example():
    for = 5  # 'for' is a reserved keyword
    print(for)

value = squareroot(16)  # squareroot doesn't exist

print(totals)  # undefined variable