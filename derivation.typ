= Model
- Pendulum consisting of $N$ parts with length 1, with masses 1 at the end. Potential contains a scale factor $g$.
- $0 <= n < N$
- Angle $phi_n$ is zero in equilibrium position (pointing towards center of gravity)
- Displacement $z_n = e^(i phi_n)$

Position of $n$th mass:

$ x_n &= sum_(j=0)^(n-1) z_j $

$
T &= 1/2 sum_n |dot(x)_n|^2
  = 1/2 sum_n sum_(j, k < n) dash(dot(z)_j) dot(z)_k
  = 1/2 sum_(j,k) (N-max(j,k)) dash(dot(z)_j) dot(z)_k
  = 1/2 sum_(j,k) (N-max(j,k)) dot(phi)_j dot(phi)_k cos(phi_j - phi_k)
\
 &= 1/2 sum_(j,k) (N-max(j,k)) c_(j k) dot(phi)_j dot(phi)_k
\
V &= -g sum_n 1/2 (x_n + dash(x)_n)
  = -g/2 sum_n (N-n) (z_n + dash(z)_n)
  = -g/2 sum_n (N-n) s_n
\
L &= T - V = 1/2 sum_(j,k) (N-max(j,k)) c_(j k) dot(phi)_j dot(phi)_k + g/2 sum_n (N-n) s_n
\
d/(d t) (partial L)/ (partial dot(phi)_n)
  &= d/(d t) ( sum_j (N-max(j,n)) c_(j n) dot(phi)_j )
   = sum_j (N-max(j,n)) (c_(j n) dot.double(phi)_j + s_(j n) (dot(phi)_n - dot(phi)_j) dot(phi)_j)
\
(partial L)/(partial phi_n)
  &= sum_j (N-max(j,n)) s_(j n) dot(phi)_j dot(phi)_n - g/2 (N-n) c_n
$
$
sum_j (N-max(j, n))c_(j n) dot.double(phi)_j
&= sum_j (N-max(j, n)) s_(j n) dot(phi)_j^2 - g/2 (N-n) c_n
$
