#include "m_pd.h"
#include "math.h"

static t_class *duck_tilde_class;

typedef struct _duck_tilde {
  t_object  x_obj;

  t_sample f;

  float power;
  float amplitude;
  float time;

  t_inlet *x_in2;
  t_inlet *x_in3;
  t_inlet *x_in4;
  t_outlet*x_out;


} t_duck_tilde;

t_int *duck_tilde_perform(t_int *w)
{
  t_duck_tilde *x = (t_duck_tilde *)(w[1]);
  t_sample  *in1 =    (t_sample *)(w[2]);
  t_sample  *in2 =    (t_sample *)(w[3]);
  t_sample  *out =    (t_sample *)(w[4]);
  int          n =           (int)(w[5]);

  x->power = x->power * x->time;
  int i;
  for (i=0;i<n;i++)
  {
      x->power += (1-x->time)*pow(*in1++,2)/n;
  }

  while (n--) *out++ = *in2++ * 1/(1+ x->amplitude*x->power);

  return (w+6);
}

void duck_tilde_dsp(t_duck_tilde *x, t_signal **sp)
{
  dsp_add(duck_tilde_perform, 5, x,
          sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
}

void duck_tilde_free(t_duck_tilde *x)
{
  inlet_free(x->x_in2);
  outlet_free(x->x_out);
}

void *duck_tilde_new(t_floatarg f)
{
  t_duck_tilde *x = (t_duck_tilde *)pd_new(duck_tilde_class);

  x->x_in2=inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  x->x_in3=floatinlet_new(&x->x_obj, &x->amplitude);
  x->x_in4=floatinlet_new(&x->x_obj, &x->time);
  x->x_out=outlet_new(&x->x_obj, &s_signal);

  return (void *)x;
}

void duck_tilde_setup(void) {
  duck_tilde_class = class_new(gensym("duck~"),
        (t_newmethod)duck_tilde_new,
        0, sizeof(t_duck_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT, 0);

  class_addmethod(duck_tilde_class,
        (t_method)duck_tilde_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(duck_tilde_class, t_duck_tilde, f);
}
