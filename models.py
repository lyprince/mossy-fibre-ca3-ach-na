from brian2 import NeuronGroup,Synapses
from brian2.units import *

########################
### CA3 Neuron Model ###
########################

def ca3_neurongroup(N,ACh=0.0):
    
    if type(ACh) is list:
        ACh_ex,ACh_syn = ACh
    else:
        ACh_ex = ACh_syn = ACh
    
    # Izhikevich Model

    eqs_iz = '''
    dv/dt = (k*(v-vt)*(v-vr)-((ge_tot/(1+ACh_syn))+geff_tot + gje_tot)*(v-ve)-(gi_tot+giff_tot)*(v-vi) + Iapp - u)/Cm: volt (unless refractory)
    du/dt = a*(b*(v-vr) - u) : amp (unless refractory)
    k : amp/volt**2
    vt : volt
    vr : volt
    ve : volt
    vi : volt
    vpeak : volt
    a : Hz
    b : siemens
    c : volt
    d : amp
    Cm : farad
    gje_tot : siemens
    ge_tot : siemens
    gi_tot : siemens
    geff_tot : siemens
    giff_tot : siemens
    Iapp : amp
    Isyn = -((ge_tot/(1+ACh_syn))+geff_tot + gje_tot)*(v-ve)-(gi_tot+giff_tot)*(v-vi) : amp
    ACh_ex : 1
    ACh_syn : 1
    '''

    iz_reset='''
    v=c
    u+=d
    '''
    
    # Create NeuronGroup
    
    nrns = NeuronGroup(N,eqs_iz, threshold='v>vpeak',reset=iz_reset,refractory=1*ms,method='rk4')
    
    # Set Parameters for CA3 Neuron from Hummos et al (2014)
    
    nrns.ACh_syn = ACh_syn
    nrns.ACh_ex = ACh_ex
    nrns.Cm = 24.*pF
    nrns.k = 1.5*pA/mV**2
    nrns.a = 10.*Hz
    nrns.b = 2.*nsiemens
    nrns.c = (-63. + 2.*nrns.ACh_ex)*mV
    nrns.d = (60. - 10.*nrns.ACh_ex)*pA
    nrns.vr = (-75.+5.*nrns.ACh_ex)*mV
    nrns.vt = -58.*mV
    nrns.ve = 10.*mV
    nrns.vi = -80.*mV
    nrns.vpeak = 29.*mV
    
    return nrns

#########################
### FSBC neuron model ###
#########################

def fsbc_neurongroup(N,ACh=0.0):
    
    # Izhikevich Model

    eqs_iz = '''
    dv/dt = (k*(v-vt)*(v-vr)-(ge_tot+geff_tot)*(v-ve)-gi_tot*(v-vi) + Iapp - u)/Cm: volt (unless refractory)
    du/dt = a*(b*(v-vr) - u) : amp (unless refractory)
    k : siemens/volt
    vt : volt
    vr : volt
    ve : volt
    vi : volt
    vpeak : volt
    a : Hz
    b : siemens
    c : volt
    d : amp
    Cm : farad
    ge_tot : siemens
    gi_tot : siemens
    geff_tot : siemens
    Iapp : amp
    ACh : 1
    '''
    iz_reset='''
    v=c
    u+=d
    '''
    
    # Create NeuronGroup
    
    nrns = NeuronGroup(N,eqs_iz, threshold='v>vpeak',reset=iz_reset,refractory=1*ms,method='rk4')
    
    # Set Parameters for CA3 Neuron from Hummos et al (2014)
    
    nrns.ACh = ACh
    nrns.Cm = 16.*pF
    nrns.k = 1.5*nS/mV
    nrns.a = 900.*Hz
    nrns.b = 2.*nsiemens
    nrns.c = -80.*mV
    nrns.d = 400.*pA
    nrns.vr = (-65.+2.*nrns.ACh)*mV
    nrns.vt = -50.*mV
    nrns.ve = 10.*mV
    nrns.vi = -80.*mV
    nrns.vpeak  = 28.*mV
    
    return nrns

###################################
### Plastic Excitatory Synapses ###
###################################

def ePlas_synapses(pre_ng,post_ng):

    esyn= '''
    dg/dt = -g/(10.*ms) : siemens (clock-driven)
    dApre/dt = -Apre/(20.*ms) : 1 (event-driven)
    dApost/dt = -Apost/(20.*ms) : 1 (event-driven)
    dAslow/dt = -Aslow/(100.*ms) : 1 (event-driven)
    dB/dt = -B/(1.*second) : 1 (event-driven)
    ge_tot_post = g : siemens (summed)
    w : 1
    wmax : 1
    w_target : 1
    A : 1
    post_rate : Hz
    '''

    esyn_pre={
        'pre_transmission' : '''g+=w*nsiemens''',
        'pre_plasticity': '''
        Apre+=1
        w=clip(w+A*(Apost*Aslow-Aslow*B**2),0.0,wmax)
        '''
    }

    esyn_post={
        'post_plasticity' : '''
        Apost+=1
        Aslow+=1
        w=clip(w+A*(Apre*Aslow-Aslow*B**2),0.0,wmax)
        B+=1./(post_rate*1.*second)
        '''
    }
    
    syns=Synapses(pre_ng,post_ng,model=esyn,on_pre=esyn_pre,on_post=esyn_post,method='linear')
    
    syns.pre_transmission.delay = 2*ms
    syns.pre_transmission.order=-2
    syns.pre_plasticity.order = -1
    syns.post_plasticity.order = -1
    
    return syns

###################################
### Plastic Inhibitory Synapses ###
###################################

def iPlas_synapses(pre_ng,post_ng):

    isyn= '''
    dg/dt = -g/(10.*ms) : siemens (clock-driven)
    dApre/dt = -Apre/(20*ms) : 1 (event-driven)
    dApost/dt = -Apost/(20*ms) : 1 (event-driven)
    gi_tot_post = g : siemens (summed)
    w : 1
    wmax : 1
    A : 1
    B : 1
    '''

    isyn_pre={
        'pre_transmission' : '''g+=w*nsiemens''',
        'pre_plasticity': '''
        Apre+=1
        w=clip(w+A*(Apost-B),0.0,wmax)
        '''
    }

    isyn_post={
        'post_plasticity' : '''
        Apost+=1
        w=clip(w+A*(Apre-B),0.0,wmax)
        '''
    }
    
    syns = Synapses(pre_ng,post_ng,model=isyn,on_pre=isyn_pre,on_post=isyn_post,method='linear')
    syns.pre_transmission.delay = 2*ms
    syns.pre_transmission.order= -2
    syns.pre_plasticity.order = -1
    syns.post_plasticity.order = -1
    
    return syns

##################################
### Static Inhibitory Synapses ###
##################################

def eStat_synapses(pre_ng,post_ng):
    
    esyn='''
    dg/dt = -g/(10.*ms) : siemens (clock-driven)
    ge_tot_post = g : siemens (summed)
    w: 1
    '''
    
    esyn_pre='''
    g+=w*nsiemens
    '''
    
    syns = Synapses(pre_ng,post_ng,model=esyn,on_pre=esyn_pre,method='linear',delay=2*ms)
    return syns

##################################
### Static Inhibitory Synapses ###
##################################

def iStat_synapses(pre_ng,post_ng):
    
    isyn='''
    dg/dt = -g/(10.*ms) : siemens (clock-driven)
    gi_tot_post = g : siemens (summed)
    w: 1
    '''
    
    isyn_pre='''
    g+=w*nsiemens
    '''
    
    syns = Synapses(pre_ng,post_ng,model=isyn,on_pre=isyn_pre,method='linear',delay=2*ms)
    return syns


##############################
### Mossy fibre Excitation ###
##############################

def eMF_synapses(source,post_ng):
    
    mf = '''
    dg/dt = -g/(10.*ms) : siemens (clock-driven)
    dp/dt = (p0-p)/taup : 1 (event-driven)
    geff_tot_post=g: siemens (summed)
    taup = 3.3*second : second
    p0 : 1
    w : 1
    '''
    
    mf_pre='''
    g+=w*(p**2)*nsiemens
    p+=0.15*(1-p)
    '''
    
    syns = Synapses(source,post_ng,model=mf,on_pre=mf_pre,method='linear')
    return syns


##############################
### Mossy Fibre Inhibition ###
##############################

def iMF_synapses(source,post_ng):
    
    mf='''
    dg/dt = -g/(20.*ms): siemens (clock-driven)
    dp/dt = (p0-p)/(1.4*second) : 1 (event-driven)
    dq/dt = (1.0-q)/(0.8*second) : 1 (clock-driven)
    dy/dt = (y0-y)/(8.0*second) : 1 (event-driven)
    giff_tot_post=g : siemens (summed)
    w: 1
    p0 : 1
    y0 : 1
    '''
    mf_pre='''
    g+=w*p*q*nsiemens
    q-=p*q
    p+=y*(1-p)
    y+=0.11*(1-y)
    '''
    
    syns=Synapses(source,post_ng,model=mf,on_pre=mf_pre,method='linear',delay=5*ms)
    return syns