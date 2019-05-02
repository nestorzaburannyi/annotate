#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <errno.h>
#include <sys/select.h>
#include <string.h>
#include <mpi.h>

//#define SERVER_PORT  12345
int SERVER_PORT=12345;

#define TRUE             1
#define FALSE            0
#define WRKTAG           0
#define IDLTAG           3

int microtimer_(int *seconds, int *microseconds)
{
        struct timeval timeout;
        int rc;
        timeout.tv_sec  = *seconds;
        timeout.tv_usec = *microseconds;
        /* use select() as timer to sleep seconds/microseconds */
        rc=select(1,NULL,NULL,NULL,&timeout);
        return(rc);
}

int match(char *regexp, char *text)
{
    if (regexp[0] == '^')
        return matchhere(regexp+1, text);
    do {    /* must look even if string is empty */
        if (matchhere(regexp, text))
            return 1;
    } while (*text++ != '\0');
    return 0;
}


/* matchhere: search for regexp at beginning of text */
int matchhere(char *regexp, char *text)
{
    if (regexp[0] == '\0')
        return 1;
    if (regexp[1] == '*')
        return matchstar(regexp[0], regexp+2, text);
    if (regexp[0] == '$' && regexp[1] == '\0')
        return *text == '\0';
    if (*text!='\0' && (regexp[0]=='.' || regexp[0]==*text))
        return matchhere(regexp+1, text+1);
    return 0;
}

/* matchstar: search for c*regexp at beginning of text */
int matchstar(int c, char *regexp, char *text)
{
    do {    /* a * matches zero or more instances */
        if (matchhere(regexp, text))
            return 1;
    } while (*text != '\0' && (*text++ == c || c == '.'));
    return 0;
}

int socket_loop_(int *x)
{
   int    i, len, rc, on = 1, master=1, irank, nproc;
   int    max_sd, new_sd;
   int    desc_ready, end_server = FALSE, shutdown_received = FALSE;
   int    close_conn;
   const int MAXPROC=1024; /* maximum number of MPI processes */
   const int BUFSIZE=64000; /* max length of request */
   const int BUFSIZE2=128000; /* split message fits two buffers */
   const int BUFSIZE3=64000000; /* output message size */
   char   *buffer,qseq[BUFSIZE2+1];
   struct timeval       timeout;
   fd_set master_set;
   fd_set working_set;
   char *pch;
   static const char shutdown_msg[] = "SHUTDOWN!";
   const int MAX_MAXSD=128;
   char m[MAX_MAXSD][BUFSIZE+1];
   int complete_flag,l;
   int    listen_sd;
   struct sockaddr_in   addr;
   int iquery,lmsg;
   char *resultbuff, *outbuff, *request; /* H=1000, MAXRES=64k; H*MAXRES */
   MPI_Status status,statuslist[MAXPROC];
   MPI_Request requestlist[MAXPROC];

   SERVER_PORT=*x;
   fprintf(stderr,"Server port set to %d from %d\n",SERVER_PORT,*x);

   /*************************************************************/
   /* Create an AF_INET stream socket to receive incoming       */
   /* connections on                                            */
   /*************************************************************/
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   fprintf(stderr,"Start socket_loop master=%i nproc=%i\n",master,nproc);

   /* initialize */
   for(i=0;i<MAX_MAXSD;i++) { memset(m[i],'\0',BUFSIZE+1); }
   resultbuff=malloc(BUFSIZE3);
   outbuff=malloc(BUFSIZE3);
   buffer=malloc(BUFSIZE+1);
   request=malloc(BUFSIZE+1);
   memset(resultbuff,0,BUFSIZE3);
   memset(outbuff,0,BUFSIZE3);
   memset(qseq,'\0',BUFSIZE2+1); 
   memset(buffer,'\0',BUFSIZE+1); 
   memset(request,'\0',BUFSIZE+1); 
   listen_sd = socket(AF_INET, SOCK_STREAM, 0);
   if (listen_sd < 0)
   {
      perror("socket() failed");
      return(-1);
   }

   /*************************************************************/
   /* Allow socket descriptor to be reuseable                   */
   /*************************************************************/
   rc = setsockopt(listen_sd, SOL_SOCKET,  SO_REUSEADDR,
                   (char *)&on, sizeof(on));
   if (rc < 0)
   {
      perror("setsockopt() failed");
      close(listen_sd);
      return(-1);
   }

   /*************************************************************/
   /* Set socket to be non-blocking.  All of the sockets for    */
   /* the incoming connections will also be non-blocking since  */
   /* they will inherit that state from the listening socket.   */
   /*************************************************************/
   rc = ioctl(listen_sd, FIONBIO, (char *)&on);
   if (rc < 0)
   {
      perror("ioctl() failed");
      close(listen_sd);
      return(-1);
   }

   /*************************************************************/
   /* Bind the socket                                           */
   /*************************************************************/
   memset(&addr, 0, sizeof(addr));
   addr.sin_family      = AF_INET;
   addr.sin_addr.s_addr = htonl(INADDR_ANY);
   addr.sin_port        = htons(SERVER_PORT);
   fprintf(stderr,"Server port = %i\n",SERVER_PORT);
   rc = bind(listen_sd,
             (struct sockaddr *)&addr, sizeof(addr));
   if (rc < 0)
   {
      perror("bind() failed");
      close(listen_sd);
      exit(-1); // want to crash server
      return(-1);
   }

   /*************************************************************/
   /* Set the listen back log                                   */
   /*************************************************************/
   rc = listen(listen_sd, 32);
   if (rc < 0)
   {
      perror("listen() failed");
      close(listen_sd);
      return(-1);
   }

   /*************************************************************/
   /* Initialize the master fd_set                              */
   /*************************************************************/
   fprintf(stderr,"Initialize socket_loop\n");
   FD_ZERO(&master_set);
   max_sd = listen_sd;
   FD_SET(listen_sd, &master_set);
   iquery=0;

   /*************************************************************/
   /* Loop waiting for incoming connects or for incoming data   */
   /* on any of the connected sockets.                          */
   /*************************************************************/
   do
   {
        if(shutdown_received == TRUE) end_server = TRUE;

      /**********************************************************/
      /* Copy the master fd_set over to the working fd_set.     */
      /**********************************************************/
      memcpy(&working_set, &master_set, sizeof(master_set));

      /**********************************************************/
      /* Rewind clock to 5 minutes.  If shutdown flag is up and */
      /* no activity after 5 minutes this program will end.     */
      /**********************************************************/
      timeout.tv_sec  = 5 * 60;
      timeout.tv_usec = 0;
      /**********************************************************/
      /* Call select() and wait 5 minutes for it to complete.   */
      /**********************************************************/
      rc = select(max_sd + 1, &working_set, NULL, NULL, &timeout);

      /**********************************************************/
      /* Check to see if the select call failed.                */
      /**********************************************************/
      if (rc < 0)
      {
         perror("  select() failed");
         break;
      }

      /**********************************************************/
      /* Check to see if the 5 minute time out expired.         */
      /**********************************************************/
      if (rc == 0 && shutdown_received == TRUE)
      {
         fprintf(stderr,"  select() timed out.  End program.\n");
         break;
      }

      /**********************************************************/
      /* One or more descriptors are readable.  Need to         */
      /* determine which ones they are.                         */
      /**********************************************************/
      desc_ready = rc;
      for (i=0; i <= max_sd  &&  desc_ready > 0; ++i)
      {
         /*******************************************************/
         /* Check to see if this descriptor is ready            */
         /*******************************************************/
         if (FD_ISSET(i, &working_set))
         {
            /****************************************************/
            /* A descriptor was found that was readable - one   */
            /* less has to be looked for.  This is being done   */
            /* so that we can stop looking at the working set   */
            /* once we have found all of the descriptors that   */
            /* were ready.                                      */
            /****************************************************/
            desc_ready -= 1;

            /****************************************************/
            /* Check to see if this is the listening socket     */
            /****************************************************/
            if (i == listen_sd)
            {
               //fprintf(stderr,"  Listening socket is readable\n");
               /*************************************************/
               /* Accept all incoming connections that are      */
               /* queued up on the listening socket before we   */
               /* loop back and call select again.              */
               /*************************************************/
               do
               {
                  /**********************************************/
                  /* Accept each incoming connection.  If       */
                  /* accept fails with EWOULDBLOCK, then we     */
                  /* have accepted all of them.  Any other      */
                  /* failure on accept will cause us to end the */
                  /* server.                                    */
                  /**********************************************/
                  new_sd = accept(listen_sd, NULL, NULL);
                  if (new_sd < 0)
                  {
                     if (errno != EWOULDBLOCK)
                     {
                        perror("  accept() failed");
                        end_server = TRUE;
                     }
                     break;
                  }
                  /**********************************************/
                  /* Add the new incoming connection to the     */
                  /* master read set                            */
                  /**********************************************/
                  fprintf(stderr,"  New incoming connection - %d\n", new_sd);
                  FD_SET(new_sd, &master_set);
                  if (new_sd > max_sd)
                     max_sd = new_sd;

                  /**********************************************/
                  /* Loop back up and accept another incoming   */
                  /* connection                                 */
                  /**********************************************/
               } while (new_sd != -1);
            }

            /****************************************************/
            /* This is not the listening socket, therefore an   */
            /* existing connection must be readable             */
            /****************************************************/
            else
            {
               //fprintf(stderr,"  Descriptor %d is readable\n", i);
               close_conn = FALSE;
               /*************************************************/
               /* Receive all incoming data on this socket      */
               /* before we loop back and call select again.    */
               /*************************************************/
               do
               {
                  /**********************************************/
                  /* Receive data on this connection until the  */
                  /* recv fails with EWOULDBLOCK.  If any other */
                  /* failure occurs, we will close the          */
                  /* connection.                                */
                  /**********************************************/
                  memset(buffer,'\0',BUFSIZE);
                  //fprintf(stderr,"#socket: receive next\n");
                  rc = recv(i, buffer, BUFSIZE, 0);
                  //fprintf(stderr,"#socket:  %d bytes received in buffer\n", rc);
                  if (rc < 0)
                  {
                     if (errno != EWOULDBLOCK)
                     {
                        perror("  recv() failed");
                        close_conn = TRUE;
                     }
                     break;
                  }

                  /**********************************************/
                  /* Check to see if the connection has been    */
                  /* closed by the client                       */
                  /**********************************************/
                  if (rc == 0)
                  {
                     fprintf(stderr,"  Connection closed\n");
                     close_conn = TRUE;
                     break;
                  }

                  /**********************************************/
                  /* Data was recevied                          */
                  /**********************************************/
                  //fprintf(stderr,"#socket:  %d bytes received in buffer\n", rc);
                  /**********************************************/
                  /* prepend last token from previous    buffer */
                  /**********************************************/
		  memset(qseq,'\0',BUFSIZE2); /* buffer must be large enough */
                  l=strlen(m[i]);
                  strncpy(qseq,m[i],l); // init qseq with previous incomplete request
                  memmove(qseq+l,buffer,strlen(buffer)); /* append buffer to qseq */
                  /**********************************************/
                  /* extract all complete requests from buffer  */
                  /**********************************************/
                  pch = strtok (qseq,"\n"); /* pch is pointer to char array */
                  while (pch != NULL)
                  {

                    complete_flag=match("</QUERY>",pch);

                    //fprintf(stderr,"do|%i flag=%i>>%s<<\n",i,complete_flag,pch);

                    /* handle exceptions */
                    if(memcmp(shutdown_msg,pch,9) == 0)
                    {
                        shutdown_received = TRUE; /* flag, no action yet! */
                        fprintf(stderr,"shutdown received\n");
                    }
                    memset(request,'\0',BUFSIZE);
                    lmsg=strlen(pch);
                    strncpy(request,pch,lmsg);
                    if(complete_flag)
                {

                    /* send complete request to all nodes in farm */
                    iquery++;
                    //fprintf(stderr,"#socket: complete query qid=%i,nproc=%i,lmsg=%i\n",iquery,nproc,lmsg);

                    for (irank=1; irank<nproc; irank++)
                    {
                        //fprintf(stderr,"#socket: query %i being sent to rank %i of %i (%i chars)\n",iquery,irank,nproc,lmsg);
                        rc=MPI_Isend(&lmsg,1,MPI_INT,irank,0,MPI_COMM_WORLD,&requestlist[irank-1]);
                    }
                    /* send qseq to master only */
                    MPI_Send(request,lmsg,MPI_CHAR,1,0,MPI_COMM_WORLD);
                    // check that ISEND finished
                    MPI_Waitall(nproc-1, requestlist, statuslist);

                    /* receive result from master */
                    rc=MPI_Recv(&lmsg,1,MPI_INT,master,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                    //fprintf(stderr,"socket: received result (rc=%i) from master %i bytes\n",rc,lmsg);
                    rc=MPI_Recv(resultbuff,lmsg,MPI_CHAR,master,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                    //fprintf(stderr,"socket: received resultbuff (rc=%i) from master %i bytes\n",rc,lmsg);
                    /* send response to client */
                    memset(outbuff,'\0',sizeof(outbuff));
                    //fprintf(stderr,outbuff,"<QUERY uid=%d>\n",iquery);

                    strncat(outbuff,resultbuff,lmsg);
                    //strncat(outbuff,"\n</QUERY>\n",10);
                    lmsg=strlen(outbuff)+1; // include string terminator
                    rc=send(i,outbuff,lmsg,MSG_NOSIGNAL);
                    //fprintf(stderr,"# Server sent %i bytes result (rc=%i)\n",lmsg,rc);
                    if(rc<0)
                    {
                        if(errno != EWOULDBLOCK)
                        {
                                perror("send() failed");
                                close_conn=TRUE;
                        }
                        break;
                    }
                    if(rc==0)
                    {
                        fprintf(stderr,"  Connection closed\n");
                        close_conn=TRUE;
                        break;
                    }
                    // request was processed: initialize m[i]
                    memset(m[i],'\0',BUFSIZE);
                }  else {
                        //fprintf(stderr,"  incomplete request (%i bytes) was not sent\n",strlen(request));
                        memset(m[i],'\0',BUFSIZE); /* max message length */
                        memmove(m[i],request,strlen(request));
                }
                    // get next token (input split on newline)
                    pch = strtok (NULL, "\n");
                  }

                  if(i < MAX_MAXSD && shutdown_received == FALSE) break;
               } while (TRUE);
             /*************************************************/
               /* If the close_conn flag was turned on, we need */
               /* to clean up this active connection.  This     */
               /* clean up process includes removing the        */
               /* descriptor from the master set and            */
               /* determining the new maximum descriptor value  */
               /* based on the bits that are still turned on in */
               /* the master set.                               */
               /*************************************************/
               if (close_conn)
               {
                  close(i);
                  FD_CLR(i, &master_set);
                  if (i == max_sd)
                  {
                     while (FD_ISSET(max_sd, &master_set) == FALSE)
                        max_sd -= 1;
                  }
               }
            } /* End of existing connection is readable */
         } /* End of if (FD_ISSET(i, &working_set)) */
      } /* End of loop through selectable descriptors */

   } while (end_server == FALSE);

   /*************************************************************/
   /* Cleanup all of the sockets that are open                  */
   /*************************************************************/
   for (i=0; i <= max_sd; ++i)
   {
      if (FD_ISSET(i, &master_set))
         close(i);
   }
   fprintf(stderr,"# Exiting socket_loop\n");
}




