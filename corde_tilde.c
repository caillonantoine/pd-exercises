#include "m_pd.h"
#include <stdlib.h>
#include <math.h>

static t_class *corde_tilde_class;

typedef struct _corde_tilde{
    t_object x_obj;
    
    float nu;
    float rho;
    float position;
    float size;
    float B;
    float r;
    
    float *zeta;
    float *amp;
    
    unsigned int m;
    unsigned int cs;
    unsigned int active;
    
    t_inlet *x_in2;
    t_inlet *x_in3;
    t_inlet *x_in4;
    t_inlet *x_in5;
    t_inlet *x_in6;
    
    t_outlet *x_out;
    
} t_corde_tilde;

void corde_tilde_bang(t_corde_tilde *x){
    int i,j;
    int pos = (int)(x->position*100);
    //PROJECTION SUR LA BASE MODALE
    for (i=0; i<x->m; i++)
    {
        x->amp[i] = 0;
        for (j=0; j<x->size; j++)
        {
            x->amp[i] += sin(3.14157*(i+1)*(pos+j)/100);
            x->amp[i] += sin(3.14157*(i+1)*(pos-j)/100);
        }
        x->amp[i] /= x->size/2;
    }
    
    //AMORTISSEMENTS
    float L = .4;
    float f0 = x->nu*44100;
    float rho = 1e3*(3.14157*pow(x->r,2));
    float T0 = pow(2*L*f0,2)*rho;
    float nf = 3e-3;
    float na = .9;
    float nb = 2.5e-2;
    float pn2;
    
    for (i=0; i<x->m; i++)
    {float
        pn2 = pow((2*i-1)*3.14157/(2*L),2);
        x->zeta[i] = .5*(T0*(nf+na/(2*3.14157*f0*i))+nb*x->B*pn2)/(T0 + x->B*pn2);
    }
}

void *corde_tilde_new(t_floatarg nu){
    t_corde_tilde *x = (t_corde_tilde *)pd_new(corde_tilde_class);
    
    x->nu          = nu;
    x->m           = 60;
    x->cs          = 0 ;
    x->active      = 0;
    
    x->B           = 4e-5;
    
    x->position    = .5;
    
    x->zeta        = (float*) malloc(x->m * sizeof(float));
    x->amp         = (float*) malloc(x->m * sizeof(float));
    
    x->x_in2 = floatinlet_new(&x->x_obj, &x->nu);
    x->x_in3 = floatinlet_new(&x->x_obj, &x->B);
    x->x_in4 = floatinlet_new(&x->x_obj, &x->position);
    x->x_in5 = floatinlet_new(&x->x_obj, &x->size);
    x->x_in6 = floatinlet_new(&x->x_obj, &x->r);
    
    x->x_out = outlet_new(&x->x_obj,&s_signal);
    return (void *)x;
}

t_int *corde_tilde_perform(t_int* w){
    
    t_corde_tilde * x = (t_corde_tilde *) w[1];
    t_sample *out     = (t_sample *)      w[2];
    int n             = (int)             w[3];
    
    x->cs = (x->cs>100000)?0:x->cs;
    unsigned int i;
    
    while (n--){
        *out = 0;
        for (i = 0; i < x->m; i++)
        {
            *out += x->amp[i]*cos(2*3.14157*(i+1)*x->nu*x->cs)/(float)x->m*exp(-x->zeta[i]*(i+1)*x->nu*x->cs);
        }
        
        out++;
        x->cs++;
    }
    return (w+4);
}

void corde_tilde_dsp(t_corde_tilde *x, t_signal **sp){
    dsp_add(corde_tilde_perform,3,x,sp[0]->s_vec,sp[0]->s_n);
}


void corde_tilde_free(t_corde_tilde *x){
    outlet_free(x->x_out);
    free(x->zeta);
    free(x->amp);
}


void corde_tilde_setup(void){
    corde_tilde_class = class_new(gensym("corde~"),
                        (t_newmethod)corde_tilde_new,
                        (t_method)corde_tilde_free,
                        sizeof(t_corde_tilde),
                        CLASS_DEFAULT,
                        A_DEFFLOAT, 0);
    
    post("Antoine CAILLON - 2018 - All right reserved lol");
    
    class_addmethod(corde_tilde_class,
                    (t_method)corde_tilde_dsp,gensym("dsp"),A_CANT,0);
    
    class_addbang(corde_tilde_class,corde_tilde_bang);
}
