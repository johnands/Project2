from scitools.std import *

def extract():
    """Leser punkter fra 'data.dat' og lagrer i array u"""

    infile = open('data.dat', 'r')

    # lagre 1. linje i n
    n = int(infile.readline())
    rho_max = float(infile.readline())
    omega_r = float(infile.readline())

    # lese filen og lagre verdiene i liste u
    u1 = []
    u2 = []
    u3 = []
    for line in infile:
        words = line.split()
        u1.append(float(words[0]))
        u2.append(float(words[1]))
        u3.append(float(words[2]))
      
    # grensebetingelser
    min_max = 0.0
    u1.insert(0, min_max)
    u1.append(min_max)
    u2.insert(0, min_max)
    u2.append(min_max)
    u3.insert(0, min_max)
    u3.append(min_max)

    # konverterer listene til arrays
    u1 = array(u1)
    u2 = array(u2)
    u3 = array(u3)
    return u1, u2, u3, n, rho_max, omega_r

u1, u2, u3, n, rho_max, omega_r = extract()
rho = linspace(0, rho_max, n+2)
prob_dist1 = (abs(u1))**2
prob_dist2 = (abs(u2))**2
prob_dist3 = (abs(u3))**2

subplot(3,1,1)
plot(rho, prob_dist1, xlabel=r'$\rho$', ylabel=r'${|u^{(\lambda_0)}(\rho)|}^2$',
     title=r'Two electrons with Coulomb repulsive interaction, $\omega_r$ = %.2f' % omega_r, 
     legend='Ground state')
subplot(3,1,2)
plot(rho, prob_dist2, xlabel=r'$\rho$', ylabel=r'${|u^{(\lambda_1)}(\rho)|}^2$',
     legend='First excited state')
subplot(3,1,3)
plot(rho, prob_dist3, xlabel=r'$\rho$', ylabel=r'${|u^{(\lambda_2)}(\rho)|}^2$',
     legend='Second excited state')
hardcopy('fig5.eps')

raw_input('Press Return key to quit: ')






