# Connecting with SSH (Linux / Mac)

## Connecting with SSH

If you want to run graphical applications from outside the university’s network, or if you prefer a graphical desktop over the console, refer to the [Connecting with NX](#connecting-with-nx) section.

### From inside the university’s network

Add the following lines to your `~/.ssh/config` file. Replace `LOGIN` with the login name provided by the administrator:

```
Host esmeralda
        Hostname 10.195.17.215
        User LOGIN
        IdentityFile ~/.ssh/esmeralda
```

**Explanation:**
*   `Host esmeralda`: This line defines an alias (`esmeralda`). The options that follow this line will only be applied when `ssh` recognizes this alias.
*   `Hostname 10.195.17.215`: Specifies the IP address of the remote host.
*   `User LOGIN`: Defines the username to use for the connection.
*   `IdentityFile ~/.ssh/esmeralda`: Tells `ssh` to use the specified private key for authentication.

With this configuration, you can connect using:

```bash
ssh esmeralda
```

**Important**: The following command will **not** work as expected because the `IdentityFile` option will not be used:

```bash
ssh LOGIN@10.195.17.215
```

The alias `esmeralda` can be changed to any name you prefer.

### From outside the university’s network

To connect from outside the university’s network, you must first request VPN access from `olivier.thibault@univ-tours.fr`. Once granted, follow the instructions provided by the administrator.

Alternatively, if you prefer not to use a graphical VPN client (like `gnome` or `network-manager`), you can:

*   Install `openvpn`.
*   Download the configuration files from the provided address.
*   Extract the archive and navigate into the directory.
*   With root privileges, run:

    ```bash
    openvpn vpn-univ-TCP-443.ovpn
    ```

Once connected via VPN, refer to the instructions in the "From inside the university’s network" section above to establish your SSH connection.

