load('richelot_aux.sage')

p = 2^8*3^5 - 1

_.<I> = GF(p)[]
K.<i> = GF(p^2, modulus=I^2+1)
E0 = EllipticCurve(K, [1, 0])

Pb = E0(0)
while (2^7)*Pb == 0:
    Pb = 3^5 * E0.random_point()
Qb = Pb
while Pb.weil_pairing(Qb, 2^8)^(2^7) == 1:
    Qb = 3^5 * E0.random_point()

Pa = E0(0)
while (3^4)*Pa == 0:
    Pa = 2^8 * E0.random_point()

Qa = Pa
while Pa.weil_pairing(Qa, 3^5)^(3^4) == 1:
    Qa = 2^8 * E0.random_point()

Sa = randint(0, 3^5-1)
R = Pa + Sa * Qa
phiA = E0.isogeny(R)
EA = phiA.codomain()
phiPb = phiA(Pb)
phiQb = phiA(Qb)
# random 3^5 basis in EA
PA,QA =  EA.gens()
PA  = 2^8*PA
QA  = 2^8*QA

# E is the domain of the f-isogeny varphi_f : E -> E0
f = 2^8 - 3^5
phi_f = E0.isogenies_prime_degree(f)[0].dual()
E = phi_f.domain()

# Algorithm 1 step 1
c = f^-1 % 2^8

# Algorithm 1 step 2
Pb2 = c*phi_f.dual()(Pb)
Qb2 = c*phi_f.dual()(Qb)

# Algorithm 1 step 3
chain, codomain = Does22ChainSplit(E, EA, 3^5*Pb2,  3^5*Qb2, phiPb,  phiQb, 8)
# Checking if the isomorphism swapped the curves in the codomain
index = -1
if codomain[0].j_invariant() == E0.j_invariant():
    index = 0
else:
    index = 1
assert index != -1

# Algorithm 1 step 4
# Map the points using the computed isogeny chain
imPAA = (E(0), PA)
for f in chain:
    imPAA = f(imPAA)

# Algorithm 1 step 5
if not (3^4 * imPAA[index]).is_zero() and (3^5 * imPAA[index]).is_zero():
    assert imPAA[index].curve().j_invariant() == E0.j_invariant()
    assert imPAA[index].curve().isogeny(imPAA[index]).codomain().j_invariant() == EA.j_invariant()
    print("A generator of ker(phiA) is %",imPAA[index])
else:
    # Algorithm 1 step 6
    # Map the points using the computed isogeny chain
    imQAA = (E(0), QA)
    for f in chain:
        imQAA = f(imQAA)
    if not (3^4 * imQAA[index]).is_zero() and (3^5 * imPAA[index]).is_zero():
        assert imQAA[index].curve().j_invariant() == E0.j_invariant()
        assert imQAA[index].curve().isogeny(imQAA[index]).codomain().j_invariant() == EA.j_invariant()
        print("A generator of ker(phiA) is %",imQAA[index])

