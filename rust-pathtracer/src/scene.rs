use crate::prelude::*;

#[allow(unused)]
pub trait Scene : Sync + Send {

    fn new() -> Self where Self: Sized;

    fn hit(&self, ray: &Ray) -> Option<State>;
}
