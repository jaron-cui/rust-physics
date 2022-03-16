use nannou::prelude::*;
use std::vec::Vec;

pub trait Panel {
    fn pos() -> (f32, f32);
    fn update();
    fn draw();
}