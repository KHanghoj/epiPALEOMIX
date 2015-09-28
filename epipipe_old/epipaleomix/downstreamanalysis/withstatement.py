import gzip
class Saved():
    def __init__(self,t):
        self.t = t
    def __enter__(self):
        if self.t.endswith('.gz'):
            self.handler = gzip.open
        else:
            self.handler = open
        self.f=self.handler(self.t)
        return self.f
    def __exit__(self, type, value, traceback):
        self.f.close()

        
t='/Users/kehanghoej/Desktop/Rsscores2000_withoutABO.txt.gz'
with Saved(t) as fin:
    print(next(fin))

# print('hellowrold')
# for l in fin:
#     print l
#     break



with gzip.open(t)  as fin1:
    print(next(fin1))

print('hellowrold')
print(next(fin1))
for l in fin1:
    print l
    break

