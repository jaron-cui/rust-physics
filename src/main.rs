use nannou::prelude::*;
use std::vec::Vec;

fn main() {
    nannou::app(model)
        .update(update)
        .simple_window(view)
        .run();
}

struct Atom {
    mass: f32,
    radius: f32,
    x: f32,
    y: f32,
    xv: f32,
    yv: f32
}

impl Atom {
    fn new(mass: f32, radius: f32) -> Atom {
        Atom {mass, radius, x: 0.0, y: 0.0, xv: 0.0, yv: 0.0}
    }

    fn at(self, x: f32, y: f32) -> Atom {
        Atom {mass: self.mass, radius: self.radius, x, y, xv: self.xv, yv: self.yv}
    }

    fn moving(self, xv: f32, yv: f32) -> Atom {
        Atom {mass: self.mass, radius: self.radius, x: self.x, y: self.y, xv, yv}
    }

    fn step(&mut self, time: f32) {
        self.x += self.xv * time;
        self.y += self.yv * time;
    }

    fn accel(&mut self, xf: f32, yf: f32) {
        self.xv += xf / self.mass;
        self.yv += yf / self.mass;
    }

    fn damp_accel(&mut self, xf: f32, yf: f32) {
        let xa = xf / self.mass;
        let ya = yf / self.mass;
        if xa.abs() > 0.01 { self.xv += xa; }
        if ya.abs() > 0.01 { self.yv += ya; }
    }

    fn collides_in(atom1: &Atom, atom2: &Atom) -> f32 {
        let dxi = atom1.x - atom2.x;
        let dxv = atom1.xv - atom2.xv;
        let dyi = atom1.y - atom2.y;
        let dyv = atom1.yv - atom2.yv;
        let d = atom1.radius + atom2.radius;

        let a = dxv * dxv + dyv * dyv;
        let b = 2.0 * (dxi * dxv + dyi * dyv);
        let c = dxi * dxi + dyi * dyi - d * d;

        let differential = b * b - 4.0 * a * c;
        if differential < 0.0 {
            return -1.0;
        }

        let sqrt = f32::sqrt(differential);
        let plus = (sqrt - b) / (2.0 * a);
        let minus = (-sqrt - b) / (2.0 * a);
        if plus < 0.0 { return minus; }
        if minus < 0.0 { return plus; }
        return if plus < minus { plus } else { minus };
    }

    fn gravity_forces(atom1: &Atom, atom2: &Atom) -> (f32, f32) {
        let dx = atom2.x - atom1.x;
        let dy = atom2.y - atom1.y;
        let distance = f32::sqrt(dx * dx + dy * dy);
        if (distance < atom1.radius + atom2.radius + 1.0) { return (0.0, 0.0); }

        let t = dx.abs() + dy.abs();
        let force = if distance.abs() < 0.5 { 0.0 } else { 40.0 * atom1.mass * atom2.mass / (distance * distance) };
        let xc = (dx / t) * force;
        let yc = (dy / t) * force;

        return (xc, yc);
    }

    /*fn normal_forces(atom1: &Atom, atom2: &Atom) -> (f32, f32) {
        let dx = atom2.x - atom1.x;
        let dy = atom2.y - atom1.y;
        let distance = f32::sqrt(dx * dx + dy * dy);
        if (distance > atom1.radius + atom2.radius + 1.0) { return (0.0, 0.0); }

        let forces = collision_forces(atom1, atom2, 0.0);
        return forces;
    }*/

    fn collision_forces(atom1: &Atom, atom2: &Atom, t: f32) -> (f32, f32) {
        if t < -0.1 || t > 1.0 { return (0.0, 0.0); }

        let cx1 = atom1.x + atom1.xv * t;
        let cy1 = atom1.y + atom1.yv * t;
        let cx2 = atom2.x + atom2.xv * t;
        let cy2 = atom2.y + atom2.yv * t;

        let dx = cx2 - cx1;
        let dy = cy2 - cy1;
        let theta = dy.atan2(dx);
        let xp = theta.cos();//(dx / (dx + dy)).abs();
        let yp = theta.sin();
        let xf = atom2.xv * atom2.mass - atom1.xv * atom1.mass;
        let yf = atom2.yv * atom2.mass - atom1.yv * atom1.mass;
        //let f = (xf * xf + yf * yf).sqrt();
        
        //let xc = xf * xp * 0.8;
        //let yc = yf * yp * 0.8;

        let uv = xf * dx + yf * dy;
        let mv2 = dx * dx + dy * dy;
        let c = uv / mv2;
        let xc = c * dx;
        let yc = c * dy;

        print!("impulse: ({}, {}), angle: {}, time: {}, sin/cos: {}/{}\n", xc, yc, theta, t, xp, yp);
        return (xc, yc);
    }

    fn apply_forces(force: &dyn Fn(&Atom, &Atom) -> (f32, f32), atoms: &mut Vec<Atom>) {
        let n = atoms.len();
    
        let mut accels = vec![(0.0, 0.0); n];
        for i in 0..(n-1) {
            let atom1 = &atoms[i];
            for j in (i+1)..n {
                let atom2 = &atoms[j];
                let impulse = force(atom1, atom2);
            
                accels[i].0 += impulse.0;
                accels[i].1 += impulse.1;
                accels[j].0 -= impulse.0;
                accels[j].1 -= impulse.1;
            }
        }

        for (atom, accel) in &mut atoms.iter_mut().zip(accels.iter()) {
            atom.accel(accel.0, accel.1);
        }
    }

    fn apply_collisions(atoms: &mut Vec<Atom>) {
        let n = atoms.len();
    
        let mut accels = vec![(0.0, 0.0); n];
        let mut steps = vec! [1.0; n];
        for i in 0..(n-1) {
            let atom1 = &atoms[i];
            for j in (i+1)..n {
                let atom2 = &atoms[j];
                let t = Atom::collides_in(atom1, atom2);
                
                if t > -0.1 && t <= 1.0 {
                    steps[i] = t;
                    steps[j] = t;
                    let impulse = Atom::collision_forces(atom1, atom2, t);
                    accels[i].0 += impulse.0;
                    accels[i].1 += impulse.1;
                    accels[j].0 -= impulse.0;
                    accels[j].1 -= impulse.1;
                }
            }
        }

        for (atom, i) in &mut atoms.iter_mut().zip(0..n) {
            atom.accel(accels[i].0, accels[i].1);
            atom.step(steps[i]);
        }
    }
}

struct Model {
    atoms: Vec<Atom>
}

fn model(_app: &App) -> Model {
    let two = vec![Atom {mass: 1.0, radius: 10.0, x: 0.0, y: 0.0, xv: 0.0, yv: 0.0},
    Atom {mass: 3.0, radius: 10.0, x: 100.0, y: 0.0, xv: -1.0, yv: 0.0}];

    let system = vec![Atom {mass: 100.0, radius: 20.0, x: 0.0, y: 0.0, xv: 0.0, yv: -0.12},
    Atom {mass: 1.0, radius: 3.0, x: 100.0, y: 0.0, xv: 0.0, yv: 6.0},
    Atom {mass: 1.6, radius: 5.0, x: 300.0, y: 0.0, xv: 0.0, yv: 3.0},
    Atom {mass: 0.3, radius: 2.0, x: 320.0, y: 0.0, xv: 0.0, yv: 4.8}];

    let pong_planet = vec![Atom::new(150.0, 150.0), Atom::new(0.2, 3.0).at(0.0, 300.0),
        Atom::new(1.0, 10.0).at(0.0, 340.0).moving(0.0, 1.0)];

    let moon_planet = vec![Atom::new(150.0, 150.0), Atom::new(0.2, 3.0).at(0.0, 300.0),
        Atom::new(1.0, 10.0).at(340.0, 0.0).moving(0.0, 3.7), Atom::new(0.8, 8.0).at(-300.0, 0.0).moving(0.0, 1.0)];
/*
    let big_system = vec![Atom::new(5000.0, 60.0), Atom::new(100.0, 10.0).at(-300.0, 0.0).moving(0.0, 24.0),
        Atom::new(10.0, 4.0).at(-253.0, 0.0).moving(0.0, 32.8),];*/
        //Atom::new(10.0, 2.0).at(-410.0, 0.0).moving(0.0, 0.5)];
    let binary = vec![Atom::new(100.0, 30.0).at(-50.0, 0.0).moving(0.0, 4.0),
        Atom::new(100.0, 30.0).at(50.0, 0.0).moving(0.0, -4.0),
        Atom::new(1.0, 8.0).at(0.0, 200.0).moving(5.8, 0.0),
        Atom::new(0.5, 14.0).at(0.0, 400.0).moving(4.0, 0.0)];

    Model {
        /*atoms: vec![Atom {mass: 1.0, radius: 10.0, x: 0.0, y: 0.0, xv: 0.0, yv: 0.0},
        Atom {mass: 3.0, radius: 10.0, x: 100.0, y: 0.0, xv: -1.0, yv: 0.0}] */
        
        atoms: binary
    }
}

fn update(_app: &App, model: &mut Model, _update: Update) {
    Atom::apply_forces(&Atom::gravity_forces, &mut model.atoms);
    Atom::apply_collisions(&mut model.atoms);
}

fn view(_app: &App, _model: &Model, frame: Frame) {
    // Prepare to draw.
    let draw = _app.draw();

    // Clear the background to purple.
    draw.background().color(PLUM);

    for atom in &_model.atoms {
        draw.ellipse().color(GREEN).x_y(atom.x, atom.y).radius(atom.radius);
    }
    //draw.ellipse().color(STEELBLUE).x_y(0.0,0.0);
    //draw.text("hey").x_y(100.0, 100.0);

    draw.to_frame(_app, &frame).unwrap();
}

fn key_pressed(app: &App, model: &mut Model, key: Key) {
    
}