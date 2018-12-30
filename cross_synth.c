#include "m_pd.h"
#include "math.h"
#include "kiss_fft/kiss_fft.h"
#include "stdlib.h"

#define BUFFER_SIZE 4096

static t_class *cross_synth_tilde_class;

typedef struct _cross_synth_tilde {
  t_object  x_obj;
  t_sample f;

  t_inlet *x_in2;
  t_outlet*x_out;

  //FFT STUFF

  kiss_fft_cfg fft_cfg;
  kiss_fft_cfg ifft_cfg;

  kiss_fft_cpx  *buffer_1, *buffer_2;

  kiss_fft_cpx  *fftout_1, *fftout_2;

  kiss_fft_cpx  *ifftout;


} t_cross_synth_tilde;

t_int *cross_synth_tilde_perform(t_int *w)
{
  t_cross_synth_tilde *x = (t_cross_synth_tilde *)(w[1]);
  t_sample  *in1 =    (t_sample *)(w[2]);
  t_sample  *in2 =    (t_sample *)(w[3]);
  t_sample  *out =    (t_sample *)(w[4]);
  int          n =           (int)(w[5]);

  //BUFFER STUFF

  for (int i=0; i<BUFFER_SIZE-n; i++)
  {
    x->buffer_1[i].r = x->buffer_1[i+n].r;
    x->buffer_2[i].r = x->buffer_2[i+n].r;

  }

  for (int i=0; i<n; i++)
  {
    x->buffer_1[i+BUFFER_SIZE-n].r = in1[i];
    x->buffer_2[i+BUFFER_SIZE-n].r = in2[i];
  }

  kiss_fft(x->fft_cfg, x->buffer_1, x->fftout_1);
  kiss_fft(x->fft_cfg, x->buffer_2, x->fftout_2);

  //SPECTRAL MODIFICATIONS

  float mag;

  for (int i=0; i<BUFFER_SIZE; i++)
  {
    mag = (sqrt(pow(x->fftout_2[i].r,2) + pow(x->fftout_2[i].i,2)))/BUFFER_SIZE;
    x->fftout_1[i].r *= mag;
    x->fftout_1[i].i *= mag;
  }

  //END OF SPECTRAL MODIFICATIONS


  kiss_fft(x->ifft_cfg, x->fftout_1, x->ifftout);

  for (int i=0; i<n; i++)
  {
    out[i] = x->ifftout[BUFFER_SIZE-n+i].r/BUFFER_SIZE;
  }

  return (w+6);
}

void cross_synth_tilde_dsp(t_cross_synth_tilde *x, t_signal **sp)
{
  dsp_add(cross_synth_tilde_perform, 5, x,
          sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
}

void cross_synth_tilde_free(t_cross_synth_tilde *x)
{
  inlet_free(x->x_in2);
  outlet_free(x->x_out);
  free(x->fft_cfg);
  free(x->ifft_cfg);
  free(x->buffer_1);
  free(x->buffer_2);
  free(x->fftout_1);
  free(x->fftout_2);
  free(x->ifftout);
}

void *cross_synth_tilde_new(t_floatarg f)
{
  t_cross_synth_tilde *x = (t_cross_synth_tilde *)pd_new(cross_synth_tilde_class);

  x->x_in2=inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);

  x->x_out=outlet_new(&x->x_obj, &s_signal);

  //FFT STUFF

  x->fft_cfg  = kiss_fft_alloc(BUFFER_SIZE, 0 ,0, 0);
  x->ifft_cfg = kiss_fft_alloc(BUFFER_SIZE, 1, 0, 0);
  x->buffer_1 = (kiss_fft_cpx*) malloc(BUFFER_SIZE*sizeof(kiss_fft_cpx));
  x->buffer_2 = (kiss_fft_cpx*) malloc(BUFFER_SIZE*sizeof(kiss_fft_cpx));
  x->fftout_1 = (kiss_fft_cpx*) malloc(BUFFER_SIZE*sizeof(kiss_fft_cpx));
  x->fftout_2 = (kiss_fft_cpx*) malloc(BUFFER_SIZE*sizeof(kiss_fft_cpx));
  x->ifftout  = (kiss_fft_cpx*) malloc(BUFFER_SIZE*sizeof(kiss_fft_cpx));

  for (int i=0; i<BUFFER_SIZE; i++)
  {
    x->buffer_1[i].r = 0;
    x->buffer_1[i].i = 0;
    x->buffer_2[i].r = 0;
    x->buffer_2[i].i = 0;
    x->fftout_1[i].r = 0;
    x->fftout_1[i].i = 0;
    x->fftout_2[i].r = 0;
    x->fftout_2[i].i = 0;
    x->ifftout[i].r  = 0;
    x->ifftout[i].i  = 0;
  }

  return (void *)x;
}

void cross_synth_tilde_setup(void) {
  cross_synth_tilde_class = class_new(gensym("cross_synth~"),
        (t_newmethod)cross_synth_tilde_new,
        0, sizeof(t_cross_synth_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT, 0);

  class_addmethod(cross_synth_tilde_class,
        (t_method)cross_synth_tilde_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(cross_synth_tilde_class, t_cross_synth_tilde, f);
}
