use std::time::{Duration,Instant};
use rand::Rng;

use nalgebra::{SMatrix,SVector};
// use nalgebra_lapack::Cholesky;

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use sdl2::gfx::primitives::DrawRenderer;

const NDF: usize = 64;

type Vector = SVector<f64, NDF>;
type Matrix = SMatrix<f64, NDF, NDF>;
type StateVector = SMatrix<f64, NDF, 2>;

struct State {
    s: StateVector,
}

impl State {
    fn sv_from_halves(x: &Vector, v: &Vector) -> StateVector {
        let mut result = StateVector::zeros();
        for i in 0..NDF {
            result[(i, 0)] = x[i];
            result[(i, 1)] = v[i];
        }
        result
    }
    fn from_halves(x: &Vector, v: &Vector) -> State {
        State {s: Self::sv_from_halves(x, v)}
    }

    fn from_rand() -> State {
        let mut rng = rand::thread_rng();
        State::from_halves(
            &Vector::from_fn(|_, _| 8.*rng.gen::<f64>() - 0.5),
            &Vector::from_fn(|_, _| 4.*rng.gen::<f64>() - 0.5),
        )
    }

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
        cos_diff.cholesky().unwrap().solve_mut(&mut y);
        // Cholesky::new(cos_diff).unwrap().solve_mut(&mut y);

        State::sv_from_halves(
            &Vector::from_iterator(v.iter().cloned()),
            &y,
        )
    }

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

    fn step(&mut self) {
        let h = 0.001;
        let k1 = Self::deriv(self.s);
        let k2 = Self::deriv(self.s + 0.5*h*k1);
        let k3 = Self::deriv(self.s + 0.5*h*k2);
        let k4 = Self::deriv(self.s + h*k3);
        self.s += h/6.0 * (k1+2.*k2+2.*k3+k4);
    }

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
        }
        canvas.filled_circle(vx[0], vy[0], 5, Color::RGB(0,0,0))?;
        canvas.bezier(&vx, &vy, 2, Color::RGB(0,0,0))?;

        Ok(())
    }
}

fn main() -> Result<(), String> {
    let sdl_context = sdl2::init()?;
    let video = sdl_context.video()?;
    let width = 1024;
    let height = 800;
    let mut canvas = video
        .window("Pendulum", width, height)
        .position_centered()
        .build()
        .map_err(|e| e.to_string())?
        .into_canvas()
        .accelerated()
        .build()
        .map_err(|e| e.to_string())?;
    let mut event_pump = sdl_context.event_pump()?;

    let mut state = State::from_rand();
    //let mut dir = 0;


    'running: loop {
        canvas.set_draw_color(Color::RGB(155, 155, 155));
        canvas.clear();
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} |
                Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                    break 'running
                },
                /*
                Event::KeyDown { keycode: Some(Keycode::Up), .. } => {
                    dir = 1;
                },
                Event::KeyDown { keycode: Some(Keycode::Down), .. } => {
                    dir = -1;
                },
                */
                _ => {}
            }
        }

        let t = Instant::now();
        for _ in 0..100 {
            state.step();
            //dir = 0;
        }
        let t1 = t.elapsed().as_micros();
        state.draw(&mut canvas)?;
        let t2 = t.elapsed().as_micros() - t1;
        println!("{} {} {}", t1, t2, state.energy());
        canvas.present();
        std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
    };
    Ok(())
}
