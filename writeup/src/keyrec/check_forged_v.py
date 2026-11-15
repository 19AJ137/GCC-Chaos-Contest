# real-data check of the keys equivalent to v: use each as the actual key of the v layer (canonical gcc.py),
# the k layer is applied from its recovered state (../verify.py), and compare with the contest files and plaintexts
import sys, os, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'gcc6')); sys.path.insert(0, os.path.join(HERE, '..'))   # gcc.py と verify.py
import gcc, verify
cdir = os.path.join(HERE, '..', '..', 'ciphertext'); pdir = os.path.join(HERE, '..', '..', 'plaintext')
A = open(os.path.join(cdir, 'A.bin'), 'rb').read(); B = open(os.path.join(cdir, 'B.bin'), 'rb').read(); E = open(os.path.join(cdir, 'E.bin'), 'rb').read()
PA = open(os.path.join(pdir, 'A.txt'), 'rb').read(); PE = open(os.path.join(pdir, 'E.txt'), 'rb').read()
print('plain_A sha256', hashlib.sha256(PA).hexdigest()[:16], 'plain_E', hashlib.sha256(PE).hexdigest()[:16])
keys = [bytes.fromhex(l.split('|')[0].strip()) for l in open(os.path.join(HERE, 'forged_v_keys.txt'), encoding='utf-8') if not l.startswith('#')]
K = verify.K; ok = 0; first = True
for key in keys:
    r = [gcc.file_cryption(key, verify.file_cryption(K, A, False), encrypt=False) == PA,
         gcc.file_cryption(key, verify.file_cryption(K, E, False), encrypt=False) == PE,
         verify.file_cryption(K, gcc.file_cryption(key, B, encrypt=False), False) == PA,
         verify.file_cryption(K, gcc.file_cryption(key, PA, encrypt=True), True) == A,
         verify.file_cryption(K, gcc.file_cryption(key, PE, encrypt=True), True) == E,
         gcc.file_cryption(key, verify.file_cryption(K, PA, True), encrypt=True) == B]
    ok += all(r)
    if first or not all(r): print(key.hex(), 'A->plain %s, E->plain %s, B->plain %s, plain->A %s, plain->E %s, plain->B %s' % tuple(r)); first = False
print('keys %d, all six checks pass for %d' % (len(keys), ok))
