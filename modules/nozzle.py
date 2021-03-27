import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import math

from turtle import Screen, Turtle, Vec2D
import sdxf
import os

def read_dat():
    headers = ['MFLOWR','CCTEMP','CCPRES','AMPRES',
                'EXPRES','K','M']
    params = [0,0,0,0,0,0,0,0,0,0,0,0]

    dat = open('modules/tfile.dat', "r")
    content = dat.readlines()
    dat.close()

    os.remove("modules/tfile.dat") 

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

def curvepoint(theta, Rt, Re):
    pos = []

    sq_eps = Re / Rt
    l_N = 0.8 * ((sq_eps - 1) * Rt) / (math.tan(math.radians(15)))

    while theta < -90:
        x = 1.5 * Rt * math.cos(math.radians(theta)) * 1000
        y = (1.5 * Rt * math.sin(math.radians(theta)) + 1.5 * Rt + Rt) * 1000

        pos.append([x, y])

        theta += 0.5

    theta = -90

    while theta < (33 - 90):
        x = 0.382 * Rt * math.cos(math.radians(theta)) * 1000
        y = (0.382 * Rt * math.sin(math.radians(theta)) + 0.382 * Rt + Rt) * 1000

        pos.append([x, y])

        theta += 0.5

    Nx = pos[-1][0]
    Ny = pos[-1][1]

    Ex = l_N * 1000
    Ey = Re * 1000

    Qx = ((Ey - math.tan(math.radians(7)) * Ex) - (Ny - math.tan(math.radians(33)) * Nx)) / (math.tan(math.radians(33)) - math.tan(math.radians(7)))
    Qy = (math.tan(math.radians(33)) * (Ey - math.tan(math.radians(7)) * Ex) - math.tan(math.radians(7)) * (Ny - math.tan(math.radians(33)) * Nx)) / (math.tan(math.radians(33)) - math.tan(math.radians(7)))

    t = 0

    while t < 1:
        x = ((1 - t)**2 * Nx + 2 * (1 - t) * t * Qx + t**2 * Ex)
        y = ((1 - t)**2 * Ny + 2 * (1 - t) * t * Qy + t**2 * Ey)

        pos.append([x, y])

        t += 0.01

    return pos

def build(Rt, Re, angle):

    # http://aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
    #
    # https://github.com/nycresistor/SDXF

    turtle = Turtle()
    turtle.penup()

    pos = curvepoint(angle, Rt, Re)

    start = Vec2D(pos[0][0], pos[0][1])

    for position in [start]:
        turtle.goto(position)
        turtle.dot()

    turtle.pendown()

    for i in range(len(pos)):
        p = Vec2D(pos[i][0], pos[i][1])

        turtle.setheading(turtle.towards(p))
        turtle.goto(p)

    turtle.dot()
    turtle.penup()

    uppos = pos

    # Draw the bottom half

    for i in range(len(pos)):
        pos[i][1] = -pos[i][1]

    start = Vec2D(pos[0][0], pos[0][1])

    for position in [start]:
        turtle.goto(position)
        turtle.dot()

    turtle.pendown()

    for i in range(len(pos)):
        p = Vec2D(pos[i][0], pos[i][1])

        turtle.setheading(turtle.towards(p))
        turtle.goto(p)

    turtle.dot()

    screen = Screen()
    screen.exitonclick()

    dnpos = pos

    d = sdxf.Drawing()

    d.append(sdxf.PolyLine(points=uppos))
    d.append(sdxf.PolyLine(points=dnpos))

    d.saveas('hello_world.dxf')

    return pos

params = noz_vals(read_dat())

build(0.5 * params['throat diameter'][0], 0.5 * params['exit diameter'][0], -135)

