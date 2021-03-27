import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from turtle import Screen, Turtle, Vec2D
import sdxf

def read_dat(file):
    headers = ['MFLOWR','CCTEMP','CCPRES','AMPRES',
                'EXPRES','K','M']
    params = [0,0,0,0,0,0,0,0,0,0,0,0]

    with open(file) as f:
        content = f.readlines()

    for i in range(len(content)):
        line = content[i].split()

        try:
            if line[0] == headers[0]:
                params[0] = float(line[1])
                params[1] = float(line[2])
            elif line[0] == headers[1]:
                params[2] = float(line[1])
                params[3] = float(line[2])
            elif line[0] == headers[2]:
                params[4] = float(line[1]) * 1e6
                params[5] = float(line[2]) * 1e6
            elif line[0] == headers[3]:
                params[6] = float(line[1]) * 1e6
                params[7] = float(line[2]) * 1e6
            elif line[0] == headers[4]:
                params[8] = float(line[1]) * 1e6
                params[9] = float(line[2]) * 1e6
            elif  line[0] == headers[5]:
                params[10] = float(line[1])
            elif line[0] == headers[6]:
                params[11] = float(line[1])

        except IndexError:
            pass

    return params

def noz_vals(params, typ='arr'):
    R = 8314
    
    if typ == 'arr':
        q = params[0]
        dq = params[1]
        Tc = params[2]
        dTc = params[3]
        Pc = params[4]
        dPc = params[5]
        Pa = params[6]
        dPa = params[7]
        Pe = params[8]
        dPe = params[9]
        k = params[10]
        M = params[11]
        
    # Numerical evaluation
    
    Pt = Pc * (1 + (k - 1) / 2) ** (-k / (k-1))
    Tt = Tc * (1 / (1 + (k - 1) / 2))
    
    # Evaluating uncertainties (Ugly)
    
    sPc, sTc, sk = sp.symbols('sPc sTc sk')
    ePt = sPc * (1 + (sk - 1) / 2) ** (-sk / (sk-1))
    eTt = sTc * (1 / (1 + (sk - 1) / 2))
    
    edPt = sp.sqrt((sp.diff(ePt, sPc) * dPc)**2)
    edTt = sp.sqrt((sp.diff(eTt, sTc) * dTc)**2)
    
    dPt = edPt.evalf(subs={sk: k})
    dTt = edTt.evalf(subs={sk: k})
    
    # Numerical evaluation

    At = (q / Pt) * np.sqrt((R * Tt) / (M * k))
    Dt = np.sqrt((4 * At) / (np.pi))
    
    # Evaluating uncertainties
    
    sPt, sTt, sq, sR, sM, sAt, sAe = sp.symbols('sPt sTt sq sR sM sAt sAe')
    
    eAt = (sq / sPt) * sp.sqrt((sR * sTt) / (sM * sk))
    edAt = sp.sqrt((sp.diff(eAt, sq) * dq)**2 + (sp.diff(eAt, sPt) * dPt)**2 + (sp.diff(eAt, sTt) * dTt)**2)
    
    dAt = edAt.evalf(subs={sq: q, sPt: Pt, sR: R, sTt: Tt, sM: M, sk: k})
    
    eDt = sp.sqrt((4 * sAt) / sp.pi)
    edDt = sp.sqrt((sp.diff(eDt, sAt) * dAt)**2)
    
    dDt = edDt.evalf(subs={sAt: At})
    
    # Numerical evaluation
    
    m = np.sqrt((2 / (k - 1)) * ((Pc / Pe)**((k-1)/k) - 1))
    Ae = (At / m) * ((1 + (k - 1) / 2 * m**2)/((k + 1) / 2))**((k+1)/(2*(k-1)))
    
    De = np.sqrt((4 * Ae) / (np.pi))
    
    Ve = np.sqrt((2 * k / (k - 1)) * (R * Tc / M) * (1 - (Pe / Pc)**((k-1)/k)))
    F = q * Ve + (Pe - Pa) * Ae

    # Evaluating uncertainties
    
    sPa, sm, sVe, sPe, sAe = sp.symbols('sPa sm sVe sPe sAe')
    
    em = sp.sqrt((2 / (sk - 1)) * ((sPc / sPe)**((sk-1)/sk) - 1))
    edm = sp.sqrt((sp.diff(em, sPc) * dPc)**2 + (sp.diff(em, sPe) * dPe)**2)
    dm = edm.evalf(subs={sk: k, sPc: Pc, sPe: Pe})
    
    
    eAe = (sAt / sm) * ((1 + (sk - 1) / 2 * sm**2)/((sk + 1) / 2))**((sk+1)/(2*(sk-1)))
    edAe = sp.sqrt((sp.diff(eAe, sAt) * dAt)**2 + (sp.diff(eAe, sm) * dm)**2)
    dAe = edAe.evalf(subs={sAt: At, sm: m, sk: k})
    
    eVe = sp.sqrt((2 * sk / (sk - 1)) * (sR * sTc / sM) * (1 - (sPe / sPc)**((sk-1)/sk)))
    edVe = sp.sqrt((sp.diff(eVe, sTc) * dTc)**2 + (sp.diff(eVe, sPe) * dPe)**2 + (sp.diff(eVe, sPc) * dPc)**2)
    dVe = edVe.evalf(subs={sk: k, sR: R, sTc: Tc, sM: M, sPe: Pe, sPc: Pc})
    
    eF = sq * sVe + (sPe - sPa) * sAe
    edF = sp.sqrt((sp.diff(eF, sq) * dq)**2 + (sp.diff(eF, sVe) * dVe)**2 + (sp.diff(eF, sPe) * dPe)**2 + (sp.diff(eF, sPa) * dPa)**2 + (sp.diff(eF, sAe) * dAe)**2)
    dF = edF.evalf(subs={sq: q, sVe: Ve, sPe: Pe, sPa: Pa, sAe: Ae})

    eDe = sp.sqrt((4 * sAe) / sp.pi)
    edDe = sp.sqrt((sp.diff(eDe, sAe) * dAe)**2)
    
    dDe = edDe.evalf(subs={sAe: Ae})
    
    params = {
        'throat pressure': [Pt, dPt],
        'throat temperature': [Tt, dTt],
        'throat area' : [At, dAt],
        'throat diameter' : [Dt, dDt],
        'exit mach' : [m, dm],
        'exit area' : [Ae, dAe],
        'exit diameter' : [De, dDe],
        'exit velocity' : [Ve, dVe],
        'thrust' : [F, dF]
    }
    
    return params

def curvepoint(A, B, C, t):
    return B + (1 - t)**2 * (A - B) + t**2 * (C - B)

def build():

    # http://aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
    #
    # https://github.com/nycresistor/SDXF

    p0 = [0, 50]
    p1 = [100, 190]
    p2 = [300, 150]
    p3 = [0, 0]

    b = lambda t: p1 + (1 - t)**2 * (p0 - p1) + t**2 * (p2 - p1)

    turtle = Turtle()
    turtle.penup()

    for position in [Vec2D(p3[0], p3[1]), Vec2D(p2[0], p2[1]), Vec2D(p1[0], p1[1]), Vec2D(p0[0], p0[1])]:
        turtle.goto(position)
        turtle.dot()

    turtle.pendown()

    t = 0
    pt = []

    while t <= 1.02:
        posx = curvepoint(p0[0], p1[0], p2[0], t)
        posy = curvepoint(p0[1], p1[1], p2[1], t)

        pt.append((posx, posy))

        position = Vec2D(posx, posy)

        turtle.setheading(turtle.towards(position))
        turtle.goto(position)

        t += 0.01

    screen = Screen()
    screen.exitonclick()

    return pt


    
filename = r'D:/arsen/Documents/jhu/Missiles/builder/testcases/nozzle.txt'
params = noz_vals(read_dat(filename))
print(params)

print(build())