#!/usr/sbin/dtrace -o out.disco.stacks -s

#pragma D option ustackframes=100

profile-97
/execname == "disco" && arg1/
{
    @[ustack()] = count();
}

tick-$1
{
    exit(0);
}
