# Connecting with SSH (Windows)

## Windows

Download and install [MobaXterm](https://mobaxterm.mobatek.net/download.html).

### Generating a pair of keys

Run MobaXterm, click on `Start local terminal`, and execute the following command:

![Start local terminal](images/mobaxterm_step1.png)


```bash
ssh-keygen -t rsa -f esmeralda
```

This command will prompt you for a _passphrase_ twice and create two files: `esmeralda` (the private key) and `esmeralda.pub` (the public key). These two files must be located in the `C:\Users\XXX\Documents\MobaXterm\home\` folder (replace `XXX` with your user name).

### What to do with the keys

**Send only the public key (`esmeralda.pub`)** to one of the administrators: `yann.jullian@univ-tours.fr` or `olivier.thibault@univ-tours.fr`. Keep the keys in their generated location.

Note that every host you use to connect to the cluster will require this private key.

### Connecting with SSH

If you want to run graphical applications from outside the university’s network, or if you prefer a graphical desktop, refer to the [Connecting with NX](#connecting-with-nx) section.

####  From inside the university’s network

In the MobaXterm window:

*   Click on `Session` (upper left corner), then `SSH`.
*   Fill the `Remote Host` field with `10.195.17.215`.
*   Tick the `Specify username` box and enter the login name provided by the administrator.
*   In the `Advanced SSH settings` tab, tick the `Use private key` box and set its location to your private key file.
*   Click `OK`.
*   A new session should appear in the left panel; double-click it to connect.

![ssh](images/ssh.png)

####  From outside the university’s network

To connect from outside the university’s network, you must first request VPN access from `olivier.thibault@univ-tours.fr`. Once granted, follow the instructions provided by the administrator.

Once connected via VPN, refer to the instructions in the "From inside the university’s network" section above to establish your SSH connection.

###  File transfer

Once connected with MobaXterm, you can easily transfer files by dragging and dropping them in the left panel.

