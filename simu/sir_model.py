from pylab import *
#from matplotlib import *
from scipy.interpolate import interp1d

def simu_sir(S0,I0,R0,kII,kI,kR,kS,T):
    S = S0
    I = I0
    R = R0
    N = int(S+I+R)
    vecS = [S0]
    vecI = [I0]
    vecR = [R0]
    vecT = [0]
    while(vecT[-1] < T):
    #for t in range(T*N):
        rate = kII*S+kI*S*I/N + kR*I + kS*R
        #print(S,I,R,rate)
        vecT.append(vecT[-1] + exponential(1/rate))
        U = rand();
        if ( U < (kII*S+kI*S*I/N)/rate ):
            I += 1
            S -= 1
        elif (U < (kII*S+kI*S*I/N+kR*I)/rate):
            I -= 1
            R += 1
        else:
            S += 1
            R -= 1
        vecS.append(S)
        vecI.append(I)
        vecR.append(R)
    return(vecT,vecS,vecI,vecR)

def ode_sir(S0,I0,R0,kII,kI,kR,kS,T):
    S = S0
    I = I0
    R = R0
    N = S+I+R
    h = 0.001
    vecS = [S0]
    vecI = [I0]
    vecR = [R0]
    vecT = [0]
    for t in arange(0,T,h):
        dI = (kII*S+kI*S*I/N - kR*I)*h
        dS = (-kII*S-kI*S*I/N + kS*R)*h
        dR = -dI-dS
        I += dI
        S += dS
        R += dR
        vecT.append(t+h)
        vecS.append(S)
        vecI.append(I)
        vecR.append(R)
    return(vecT,vecS,vecI,vecR)

def exact_sir(S0,I0,R0,kII,kI,kR,kS,T):
    N = int(S0+I0+R0)
    h = 0.005
    SI = array([ [0. for I in range(N+1)] for S in range(N+1)]);
    AllI = array([ [I for I in range(N+1)] for S in range(N+1)]);
    AllS = array([ [S for I in range(N+1)] for S in range(N+1)]);
    SI[int(S0)][int(I0)] = 1;
    vecS = [S0]
    vecI = [I0]
    vecR = [R0]
    vecT = [0]
    for t in arange(0,T,h):
        print('\r{0:.2f} / {1}                '.format(t,T),end='')
        #print(SI)
        dSI = [ [0. for i in range(N+1)] for j in range(N+1)];
        for S in range(0,N+1):
            for I in range(0,N+1-S):
                if (S>0 and I < N):
                    dSI[S][I]     -= h*(kII*S + kI*S*I/N)*SI[S][I];
                    dSI[S-1][I+1] += h*(kII*S + kI*S*I/N)*SI[S][I];
                if (I>0):
                    dSI[S][I]   -= h*kR*I * SI[S][I];
                    dSI[S][I-1] += h*kR*I * SI[S][I];
                if (S < N):
                    dSI[S][I]   -= h*kS*(N-S-I) * SI[S][I];
                    dSI[S+1][I] += h*kS*(N-S-I) * SI[S][I];
        for S in range(0,N+1):
            for I in range(0,N+1-S):
                SI[S][I] += dSI[S][I]
        vecS.append( sum(SI*AllS))
        vecI.append( sum(SI*AllI))
        vecR.append( NaN )
        vecT.append(t)
    return(vecT,vecS,vecI,vecR)

def norm2_exact_sir(S0,I0,R0,kII,kI,kR,kS,T):
    (t,s,i,r) = ode_sir(S0,I0,R0,kII,kI,kR,kS,T);
    N = int(S0+I0+R0)
    h = 0.001
    SI = array([ [0. for I in range(N+1)] for S in range(N+1)]);
    SI[int(S0)][int(I0)] = 1;
    vecS = [0]
    vecT = [0]
    for (k,t) in enumerate(arange(0,T,h)):
        #print(SI)
        dSI = [ [0. for i in range(N+1)] for j in range(N+1)];
        for S in range(0,N+1):
            for I in range(0,N+1-S):
                if (S>0 and I < N):
                    dSI[S][I]     -= h*(kII*S + kI*S*I/N)*SI[S][I];
                    dSI[S-1][I+1] += h*(kII*S + kI*S*I/N)*SI[S][I];
                if (I>0):
                    dSI[S][I]   -= h*kR*I * SI[S][I];
                    dSI[S][I-1] += h*kR*I * SI[S][I];
                if (S < N):
                    dSI[S][I]   -= h*kS*(N-S-I) * SI[S][I];
                    dSI[S+1][I] += h*kS*(N-S-I) * SI[S][I];
        for S in range(0,N+1):
            for I in range(0,N+1-S):
                SI[S][I] += dSI[S][I]
        AllS = array([ [(S-s[k])**2 for I in range(N+1)] for S in range(N+1)]);
        vecS.append( sum(SI*AllS))
        vecT.append(t)
    return(vecT,vecS,NaN,NaN)

def norm2_simu_sir(S0,I0,R0,kII,kI,kR,kS,T):
    (t,s,i,r) = ode_sir(S0,I0,R0,kII,kI,kR,kS,T);
    myT = linspace(0,T,20);
    f = interp1d(t,s);
    my_s = f(myT);
    my_S = zeros(len(myT))
    nb_samples = 1000
    for i in range(nb_samples):
        (dT,S,I,R) = simu_sir(S0,I0,R0,kII,kI,kR,kS,T);
        f = interp1d(dT,S)
        my_S += (f(myT)-my_s)**2/nb_samples;
    return(myT,my_S,[],[])

def average_simu_sir(S0,I0,R0,kII,kI,kR,kS,T):
    if S0 <= 100:
        nb_samples = 10000
    else:
        nb_samples = 1000
    a = mean(array([simu_sir(S0,I0,R0,kII,kI,kR,kS,T) for i in range(0,nb_samples)]),0)
    return(a[0,:],a[1,:],a[2,:],a[3,:])



for N in [10,100]:
    for kI in [4,5]:
        (T,S,I,R) = simu_sir(.4*N,.4*N,.2*N,1.,kI,1.,1.,6)
        (t,s,i,r) = ode_sir(.4*N,.4*N,.2*N,1.,kI,1.,1.,4)
        (Te,Se,Ie,Re) = exact_sir(.4*N,.4*N,.2*N,1.,kI,1.,1.,4)
        
        f = figure(N);
        f.set_size_inches(4,4); 
        clf();
        plot(T,array(S)/N)
        plot(t,array(s)/N,'g--', linewidth=3)
        plot(Te,array(Se)/N,'r:', linewidth=3)
        xlim([0,4])
        ylim([0,0.60])
        xlabel('time');
        ylabel('X_S');
        legend(('$X^N$ (one sample path)','Mean-field approx.','$E[X^N]$ (exact)'))
        f.savefig("sir_VS_mean{}_kI{}.pdf".format(N,kI),bbox_inches='tight')
        print('done N={0}, kI={1}'.format(N,kI))

STOP
quit()

myN = range(10,100,10);
myN = [100,150,200,300];
my_s  = zeros(len(myN))
try:
    p = pd.DataFrame.from_csv('SIR_N.csv')
except:
    print("Unexpected error:", sys.exc_info()[0])
    p = pd.DataFrame({'N':[],'s':[],'Se':[],'T':[],'V':[]})
for (k,N) in enumerate(myN):
    print(N)
    if not ((p['N']==N)).any():
        (t,s,i,r) = ode_sir(.4*N,.4*N,.2*N,1.,1.,1.,1.,2.1)
        (Te,Se,Ie,Re) = exact_sir(.4*N,.4*N,.2*N,1.,1.,1.,1.,2.1)
        #(vT,vS,vI,vR) = norm2_simu_sir(.4*N,.4*N,.2*N,1.,1.,1.,1.,2)
        (vTe,vSe,_,_) = norm2_exact_sir(.4*N,.4*N,.2*N,1.,1.,1.,1.,2.1)
        for T in [.1, .2, .5, .7, 1, 1.5, 2]:
            v_s = interp1d(t,s)(T)
            v_Se = interp1d(Te,Se)(T)
            v_Ve = interp1d(vTe,vSe)(T);
            p = p.append({'N':N, 's':v_s, 'Se':v_Se, 'T':T, 'V':v_Ve}, ignore_index=True)
p.to_csv('SIR_N.csv')


fM = figure(1); clf();
fV = figure(2); clf();
for T in [.5,.7,1,2]:
    t = p['T']==T;
    figure(1); plot(p['N'][t], abs(p['s'][t]-p['Se'][t]),'+-')
    figure(2); plot(p['N'][t], p['V'][t]/p['N'][t],'+-')
figure(1); ylim([0,.08]); xlabel('N'); ylabel('$N|P(S^{(N)}_1(T) = S) - x_S(t)|$')
legend(('T=0.5','T=0.7','T=1','T=2'),loc='best')
figure(2); ylim([0,.25]); xlabel('N'); ylabel('$\sqrt{N E |X_S(t)- x_S(t)|^2}$')
legend(('T=0.5','T=0.7','T=1','T=2'),loc='best')
fM.savefig('sir_distance_timeN_t2.pdf')
fV.savefig('sir_variance_timeN_t2.pdf')


for N in [10, 100, 1000]:
    (T,S,I,R) = simu_sir(.4*N,.4*N,.2*N,1.,1.,1.,1.,10)
    (t,s,i,r) = ode_sir(.4*N,.4*N,.2*N,1.,1.,1.,1.,10)
    (Ta,Sa,Ia,Ra) = average_simu_sir(.4*N,.4*N,.2*N,1.,1.,1.,1.,10)
    f = figure();  clf()
    f.set_size_inches(6,4)
    plot(T,array(S)/N)
    plot(Ta,array(Sa)/N,'r-', linewidth=2)
    plot(t,array(s)/N,'g--', linewidth=3)
    xlim([0,8])
    ylim([0,.6])
    xlabel('time')
    ylabel('$X_S$');
    legend(('$X_S$ (one simulation)','$\mathbb{E}[X_S]$ (ave. over 10000 simu)','$x_s$ (mean-field approximation)'))
    f.savefig('sir_ODE_S_{}.pdf'.format(N),bbox_inches='tight')

    f = figure();  clf()
    f.set_size_inches(7,5)
    plot(T,array(I)/N)
    plot(Ta,array(Ia)/N,'r-', linewidth=2)
    plot(t,array(i)/N,'g--', linewidth=3)
    xlim([0,8])
    xlabel('time')
    ylabel('$X_I$');
    legend(('$X_I$ (simulation)','$\mathbb{E}[X_I]$ (ave. over 10000 simu)','$x_I$ (mean-field approximation)'),loc='best')
    f.savefig('sir_ODE_I_{}.pdf'.format(N),bbox_inches='tight')
