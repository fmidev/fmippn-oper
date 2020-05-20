
# Sets the variable "BeginTime" in ISO-8601 standard (e.g. "2020-04-29T17:04:29+0000")
set_BeginTime () {
   BeginTime=`date -u -Iseconds`
   BeginStamp=`echo ${BeginTime%+*} | tr T ' '`
}

timediff() {
   Dtime=$((`date -d "$EndTime" +%s`-`date -d "$BeginTime" +%s`))
}

# This function sets the variable "Runtime" as seconds between $BeginTime 
# and the current time. It also sets the variables below.
get_Runtime() {
   EndTime=`date -u -Iseconds`
   EndStamp=`echo ${EndTime%+*} | tr T ' '`
   timediff
   Runtime=$Dtime
}
