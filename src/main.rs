use std::time::{Duration,Instant};
use std::thread;
use rand::Rng;

use crossbeam::{channel::{bounded,Sender,Receiver}};
use nalgebra::{SMatrix,SVector};

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use sdl2::gfx::primitives::DrawRenderer;

const NDF: usize = 64;

type Vector = SVector<f64, NDF>;
type Matrix = SMatrix<f64, NDF, NDF>;
// Contains a position and a velocity vector
type StateVector = SMatrix<f64, NDF, 2>;

struct State {
    s: StateVector,
    receiver: Receiver<StateVector>,
}

impl State {
    // Create state vector from two vectors
    fn sv_from_halves(x: &Vector, v: &Vector) -> StateVector {
        let mut result = StateVector::zeros();
        for i in 0..NDF {
            result[(i, 0)] = x[i];
            result[(i, 1)] = v[i];
        }
        result
    }

    // Create state from two vectors
    fn from_halves(x: &Vector, v: &Vector) -> State {
        let s = Self::sv_from_halves(x, v);
        let s_for_worker = s.clone();

        let (sender, receiver) = bounded(16);
        thread::spawn(move || {
            State::simulate(s_for_worker, sender);
        });
        State {s, receiver}
    }

    // Random initialization
    fn from_rand() -> State {
        let mut rng = rand::thread_rng();
        State::from_halves(
            &Vector::from_fn(|_, _| 8.*rng.gen::<f64>() - 0.5),
            &Vector::from_fn(|_, _| 4.*rng.gen::<f64>() - 0.5),
        )
    }

    /* Calculate the derivative as given by Euler-Lagrange.
     *
     * The ODE essentially boils down to
     * C y'' = s(y) - S y'²
     * where S and C are matrices containing sinus and cosinus of differences of angles (see below)
     * and s contains the sinus of angles (from the potential)
     */
    fn deriv(p: StateVector) -> StateVector {
        let x = p.column(0);
        let v = p.column(1);
        let mut sin = Vector::zeros();
        let mut cos = Vector::zeros();
        for i in 0..NDF {
            (sin[i], cos[i]) = x[i].sin_cos();
        }
        // s_{jk} = (N-max(j,k))*sin(phi_j-phi_k)
        // c_{jk} = (N-max(j,k))*cos(phi_j-phi_k)
        let mut sin_diff = Matrix::zeros();
        let mut cos_diff = Matrix::zeros();
        for i in 0..NDF {
            let a = (NDF - i) as f64;
            cos_diff[(i,i)] = a;
            for j in 0..i {
                sin_diff[(i,j)] = a * (sin[i]*cos[j] - cos[i]*sin[j]);
                sin_diff[(j,i)] = -sin_diff[(i,j)];
                cos_diff[(i,j)] = a * (cos[i]*cos[j] + sin[i]*sin[j]);
                cos_diff[(j,i)] = cos_diff[(i,j)];
            }
        }
        // RHS vector
        let mut y = Vector::from_fn(|i, _| {
            5. * (NDF-i) as f64 * sin[i]
        }) - sin_diff * Vector::from_iterator(v.iter().map(|x| x*x));
        // Use cholesky decomposition to invert symmetric matrix
        cos_diff.cholesky().unwrap().solve_mut(&mut y);

        State::sv_from_halves(
            &Vector::from_iterator(v.iter().cloned()),
            &y,
        )
    }

    // Thread worker
    fn simulate(mut s: StateVector, sender: Sender<StateVector>) {
        let h = 0.001;
        while sender.send(s).is_ok() {
            for _ in 0..100 {
                let k1 = Self::deriv(s);
                let k2 = Self::deriv(s + 0.5*h*k1);
                let k3 = Self::deriv(s + 0.5*h*k2);
                let k4 = Self::deriv(s + h*k3);
                s += h/6.0 * (k1+2.*k2+2.*k3+k4);
            }
        }
    }

    fn steps(&mut self) {
        self.s = self.receiver.recv().unwrap();
    }

    // Calculate current energy, which is printed after draw() to check stability
    fn energy(&self) -> f64 {
        let x = self.s.column(0);
        let mut sin = Vector::zeros();
        let mut cos = Vector::zeros();
        for i in 0..NDF {
            (sin[i], cos[i]) = x[i].sin_cos();
        }
        let mut result = 0.;
        let mut cos_diff = Matrix::zeros();
        for i in 0..NDF {
            let a = (NDF - i) as f64;
            result += a*cos[i];
            cos_diff[(i,i)] = a;
            for j in 0..i {
                cos_diff[(i,j)] = a * (cos[i]*cos[j]+sin[i]*sin[j]);
                cos_diff[(j,i)] = cos_diff[(i,j)];
            }
        };
        let v = self.s.column(1);
        5.*result + 0.5 * (v.transpose() * (cos_diff * v))[(0,0)]
    }

    // Draw rope (bezier curve over the segments)
    fn draw(&self, canvas: &mut Canvas<Window>) -> Result<(), String> {
        let mut pos = (512., 384.);
        let mut vx = [0i16; NDF+1];
        let mut vy = [0i16; NDF+1];
        vx[0] = 512;
        vy[0] = 384;
        for i in 0..NDF {
            let (s, c) = self.s[i].sin_cos();
            pos = (
                pos.0 + 500./NDF as f64 * s,
                pos.1 - 500./NDF as f64 * c,
            );
            vx[i+1] = pos.0 as i16;
            vy[i+1] = pos.1 as i16;
            //canvas.thick_line(
            //    vx[i], vy[i], vx[i+1], vy[i+1], 2, Color::RGB(0,0,0)
            //)?;
        }
        canvas.filled_circle(vx[0], vy[0], 5, Color::RGB(0,0,0))?;
        canvas.bezier(&vx, &vy, 2, Color::RGB(0,0,0))?;

        Ok(())
    }
}

fn main() -> Result<(), String> {
    let sdl_context = sdl2::init()?;
    let video = sdl_context.video()?;
    let mut canvas = video
        .window("Pendulum", 1024, 800)
        .position_centered()
        .build()
        .map_err(|e| e.to_string())?
        .into_canvas()
        .accelerated()
        .build()
        .map_err(|e| e.to_string())?;
    let mut event_pump = sdl_context.event_pump()?;

    let mut state = State::from_rand();


    'running: loop {
        let t = Instant::now();
        canvas.set_draw_color(Color::RGB(155, 155, 155));
        canvas.clear();
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} |
                Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                    break 'running
                },
                _ => {}
            }
        }

        state.steps();
        let t1 = t.elapsed().as_micros();  // simulation time (wait for result)
        state.draw(&mut canvas)?;
        canvas.present();
        let t2 = t.elapsed().as_micros() - t1;  // draw time
        let energy = state.energy();
        let total = t.elapsed().as_nanos() as u32;
        println!("{t1:6} {t2:6} {energy}");
        if total < 1_000_000_000u32 / 60 {
            // target framerate: 60
            std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60 - total));
        }
    };
    Ok(())
}
