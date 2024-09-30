use std::time::Duration;
//use std::ops::{Index,IndexMut,Mul,Add,AddAssign};
use rand::Rng;

use nalgebra::{SMatrix,SVector};

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use sdl2::gfx::primitives::DrawRenderer;

const NDF: usize = 20;

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
        // s_{jk} = (N-max(j,k))*sin(phi_j-phi_k)
        // c_{jk} = (N-max(j,k))*cos(phi_j-phi_k)
        let mut s = Matrix::zeros();
        let mut c = Matrix::zeros();
        let x = p.column(0);
        let v = p.column(1);
        for i in 0..NDF {
            let a = (NDF - i) as f64;
            c[(i,i)] = a;
            for j in 0..i {
                let (sij, cij) = (x[i]-x[j]).sin_cos();
                s[(i,j)] = a * sij;
                s[(j,i)] = -s[(i,j)];
                c[(i,j)] = a * cij;
                c[(j,i)] = c[(i,j)];
            }
        }
        // RHS vector
        let mut y = Vector::from_fn(|i, _| {
            5. * (NDF-i) as f64 * x[i].sin()
        }) - s * Vector::from_iterator(v.iter().map(|x| x*x));
        c.cholesky().unwrap().solve_mut(&mut y);

        State::sv_from_halves(
            &Vector::from_iterator(v.iter().cloned()),
            &y,
        )
    }

    /*
    fn energy(&self) -> f64 {
        0.
    }
    */

    fn step(&mut self) {
        let h = 0.001;
        let k1 = Self::deriv(self.s);
        /*
        println!("{} {}", self.s, k1);
        panic!();
        */
        let k2 = Self::deriv(self.s + 0.5*h*k1);
        let k3 = Self::deriv(self.s + 0.5*h*k2);
        let k4 = Self::deriv(self.s + h*k3);
        self.s += h/6.0 * (k1+2.*k2+2.*k3+k4);
    }

    fn draw(&self, canvas: &mut Canvas<Window>) -> Result<(), String> {
        //println!("{}", self.energy());
        let mut pos = (512., 384.);
        canvas.filled_circle(pos.0 as i16, pos.1 as i16, 5, Color::RGB(0,0,0))?;
        for i in 0..NDF {
            let (s, c) = self.s[i].sin_cos();
            let nextpos = (
                pos.0 + 500./NDF as f64 * s,
                pos.1 - 500./NDF as f64 * c,
            );
            canvas.thick_line(
                pos.0 as i16, pos.1 as i16,
                nextpos.0 as i16, nextpos.1 as i16,
                2, Color::RGB(0,0,0),
            )?;
            pos = nextpos;
            canvas.filled_circle(
                pos.0 as i16, pos.1 as i16, 5, Color::RGB(0,0,0)
            )?;
        }
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

        for _ in 0..10 {
            state.step();
            //dir = 0;
        }
        state.draw(&mut canvas)?;
        canvas.present();
        std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
    };
    Ok(())
}
