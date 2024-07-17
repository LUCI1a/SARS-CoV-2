import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
t0 = 0
t_end = 24
t_eval = np.linspace(t0, t_end, 5000)
integ0 = [10.0,0,0,0,0,0,0,0,0,0,0,0]
k_bind,d_V,k_diss,k_fuse,k_uncoat,d_endosome,k_transl,f_ORF1,d_NSP,k_tr_neg,K_NSP,d_gRNA,d_gRNA_neg,k_tr_pos,k_complex,K_N,f_N,f_SP,d_N,d_SP,n_SP,n_N,K_Vrel,k_assemb,d_N_gRNA,k_release,d_assembled=12,0.12,0.61,0.5,0.5,0.06,45360,21000,0.069,3,100,0.2,0.1,1000,0.4,5000000,1200,10000,0.023,0.044,2000,456,1000,1,0.2,8,0.06
def lorenz(t,integ,k_bind,d_V,k_diss,k_fuse,k_uncoat,d_endosome,k_transl,f_ORF1,d_NSP,k_tr_neg,K_NSP,d_gRNA,d_gRNA_neg,k_tr_pos,k_complex,K_N,f_N,f_SP,d_N,d_SP,n_SP,n_N,K_Vrel,k_assemb,d_N_gRNA,k_release,d_assembled ):
  V_free,V_bound,V_endosome,gRNA_pos,NSP,gRNA_neg,gRNA,N,SP,N_gRNA,V_assembled,V_released = integ
  theta_RdRp=NSP/(NSP+K_NSP)
  theta_complex=N/(N+K_N)
  theta_assemb=SP/(SP+K_Vrel*n_SP)
  F1=-k_bind*V_free-d_V*V_free+k_diss*V_bound
  F2=k_bind*V_free-(k_fuse+k_diss+d_V)*V_bound
  F3=k_fuse*V_bound-(k_uncoat+d_endosome)*V_endosome
  F4=k_uncoat*V_endosome-d_gRNA*gRNA_pos
  F5=k_transl*gRNA_pos/f_ORF1-d_NSP*NSP
  F6=k_tr_neg*gRNA_pos*theta_RdRp-d_gRNA_neg*gRNA_neg
  F7=k_tr_pos*gRNA_neg*theta_RdRp-(k_complex*theta_complex+d_gRNA)*gRNA
  F8=k_transl*gRNA/f_N-k_complex*n_N*theta_complex*gRNA-d_N*N
  F9=k_transl*gRNA/f_SP-k_assemb*n_SP*theta_assemb*N_gRNA-d_SP*SP
  F10=k_complex*theta_complex*gRNA-(k_assemb*theta_assemb+d_N_gRNA)*N_gRNA
  F11=k_assemb*theta_assemb*N_gRNA-(k_release+d_assembled)*V_assembled
  F12=k_release*V_assembled-d_V*V_released
  return [F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12]

sol = solve_ivp(lorenz, [t0, t_end], integ0, t_eval=t_eval, args=(k_bind,d_V,k_diss,k_fuse,k_uncoat,d_endosome,k_transl,f_ORF1,d_NSP,k_tr_neg,K_NSP,d_gRNA,d_gRNA_neg,k_tr_pos,k_complex,K_N,f_N,f_SP,d_N,d_SP,n_SP,n_N,K_Vrel,k_assemb,d_N_gRNA,k_release,d_assembled))

fig, axs = plt.subplots(4, 1, figsize=(5, 10))
# axs[0].plot(sol.t, sol.y[0], 'b', lw=2)
# axs[0].set_ylabel('Vfree')
# axs[1].plot(sol.t, sol.y[1], 'r', lw=2)
# axs[1].set_ylabel('Vbound')
# axs[2].plot(sol.t, sol.y[2], 'g', lw=2)
# axs[2].set_ylabel('Vendsome')
# axs[3].plot(sol.t, sol.y[3], 'y', lw=2)
# axs[3].set_ylabel('gRNA(+)')

# axs[0].plot(sol.t, sol.y[4], 'b', lw=2)
# axs[0].set_ylabel('NSP')
# axs[1].plot(sol.t, sol.y[5], 'r', lw=2)
# axs[1].set_ylabel('gRNA(-)')
# axs[2].plot(sol.t, sol.y[6], 'g', lw=2)
# axs[2].set_ylabel('gRNA')
# axs[3].plot(sol.t, sol.y[7], 'y', lw=2)
# axs[3].set_ylabel('N')

axs[0].plot(sol.t, sol.y[8], 'b', lw=2)
axs[0].set_ylabel('SP')
axs[1].plot(sol.t, sol.y[9], 'r', lw=2)
axs[1].set_ylabel('N-gRNA')
axs[2].plot(sol.t, sol.y[10], 'g', lw=2)
axs[2].set_ylabel('Vassembled')
axs[3].plot(sol.t, sol.y[11], 'y', lw=2)
axs[3].set_ylabel('Vreleased')
axs[3].set_xlabel('t')
plt.show()
